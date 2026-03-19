import os
import subprocess
from pathlib import Path
from typing import Tuple, Optional

# ── PDBFixer import ──────────────────────────────────────────────────────────
try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    _PDBFIXER_AVAILABLE = True
except ImportError:
    _PDBFIXER_AVAILABLE = False

# ── Meeko import (ligand + receptor prep for AutoDock/GNINA) ─────────────────
try:
    from meeko import MoleculePreparation
    from meeko import PDBQTMolecule
    _MEEKO_AVAILABLE = True
except ImportError:
    _MEEKO_AVAILABLE = False

JOBS_DIR = Path("storage/jobs")


def prepare_receptor(pdb_path: str, job_id: str) -> str:
    """
    Full receptor preparation pipeline:
      1. PDBFixer — removes water, adds missing residues, adds hydrogens at pH 7.4
      2. Exports a clean fixed .pdb
      3. Converts to .pdbqt format for GNINA using obabel (Open Babel)

    Returns the path to the final .pdbqt receptor file.
    Raises RuntimeError if preparation fails at any step.
    """
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    fixed_pdb_path = str(job_dir / "receptor_fixed.pdb")
    pdbqt_path     = str(job_dir / "receptor.pdbqt")

    # Step 1: Fix the PDB with PDBFixer
    if _PDBFIXER_AVAILABLE:
        _fix_pdb_with_pdbfixer(pdb_path, fixed_pdb_path)
    else:
        print("[ProteinPrep] PDBFixer not available — using raw PDB (less accurate).")
        _strip_waters_manually(pdb_path, fixed_pdb_path)

    # Step 2: Convert fixed .pdb → .pdbqt using Open Babel (obabel)
    # obabel is installed in the Docker image — handles atom typing for GNINA
    _convert_pdb_to_pdbqt(fixed_pdb_path, pdbqt_path)

    return pdbqt_path


def prepare_ligand(smiles: str, job_id: str, compound_index: int) -> str:
    """
    Prepares a single ligand SMILES string for GNINA docking:
      1. Generates 3D coordinates using RDKit ETKDG
      2. Converts to .pdbqt using Meeko or obabel fallback

    Returns the path to the ligand .pdbqt file.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    job_dir = JOBS_DIR / job_id / "ligands"
    job_dir.mkdir(parents=True, exist_ok=True)

    sdf_path   = str(job_dir / f"ligand_{compound_index}.sdf")
    pdbqt_path = str(job_dir / f"ligand_{compound_index}.pdbqt")

    # Generate 3D coordinates
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot generate 3D coords: invalid SMILES '{smiles}'")

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        # ETKDGv3 failed — try random coords
        result = AllChem.EmbedMolecule(mol, AllChem.EmbedParameters())
    if result == -1:
        raise RuntimeError(
            f"RDKit could not generate 3D coordinates for SMILES: '{smiles}'. "
            "The molecule may be too rigid or have unusual valence."
        )

    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    # Write SDF
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    # Convert SDF → PDBQT
    if _MEEKO_AVAILABLE:
        _convert_sdf_to_pdbqt_meeko(mol, pdbqt_path)
    else:
        _convert_sdf_to_pdbqt_obabel(sdf_path, pdbqt_path)

    return pdbqt_path


# ── Internal helpers ─────────────────────────────────────────────────────────

def _fix_pdb_with_pdbfixer(input_pdb: str, output_pdb: str) -> None:
    """
    Uses PDBFixer to:
      - Remove all water molecules (HOH, WAT)
      - Find and fill missing residues
      - Find and fill missing heavy atoms
      - Add hydrogens at pH 7.4
    """
    fixer = PDBFixer(filename=input_pdb)

    fixer.removeHeterogens(keepWater=False)     # Strip waters + non-protein HETATM
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    with open(output_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


def _strip_waters_manually(input_pdb: str, output_pdb: str) -> None:
    """
    Minimal fallback when PDBFixer is unavailable.
    Strips water records (HOH, WAT) from the PDB file line by line.
    Does NOT add missing atoms or hydrogens.
    """
    water_res = {"HOH", "WAT", "DOD"}
    lines_out = []
    with open(input_pdb, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record in ("ATOM", "HETATM"):
                res_name = line[17:20].strip()
                if res_name in water_res:
                    continue
            lines_out.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(lines_out)


def _convert_pdb_to_pdbqt(input_pdb: str, output_pdbqt: str) -> None:
    """
    Converts a prepared receptor .pdb to .pdbqt using Open Babel.
    obabel must be installed in the Docker image (apt-get install openbabel).
    The -xr flag removes non-polar hydrogens (standard for receptor prep).
    """
    cmd = [
        "obabel",
        input_pdb,
        "-O", output_pdbqt,
        "-xr",          # remove non-polar hydrogens
        "--partialcharge", "gasteiger",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

    if result.returncode != 0:
        raise RuntimeError(
            f"Open Babel receptor conversion failed:\n"
            f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
        )
    if not Path(output_pdbqt).exists():
        raise RuntimeError(
            f"Open Babel ran but did not produce output file: {output_pdbqt}"
        )


def _convert_sdf_to_pdbqt_meeko(mol, output_pdbqt: str) -> None:
    """
    Converts an RDKit Mol object to .pdbqt using Meeko.
    Meeko correctly handles flexible torsions and atom typing for AutoDock/GNINA.
    """
    preparator = MoleculePreparation()
    preparator.prepare(mol)
    pdbqt_string = preparator.write_pdbqt_string()
    with open(output_pdbqt, "w") as f:
        f.write(pdbqt_string)


def _convert_sdf_to_pdbqt_obabel(sdf_path: str, output_pdbqt: str) -> None:
    """
    Fallback ligand conversion using Open Babel when Meeko is unavailable.
    Less accurate torsion handling than Meeko but fully functional.
    """
    cmd = [
        "obabel",
        sdf_path,
        "-O", output_pdbqt,
        "--partialcharge", "gasteiger",
        "-xh",      # add polar hydrogens
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

    if result.returncode != 0:
        raise RuntimeError(
            f"Open Babel ligand conversion failed:\n"
            f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
        )

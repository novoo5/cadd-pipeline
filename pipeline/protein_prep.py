import os
import subprocess
from pathlib import Path
from typing import Tuple, Optional

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    _PDBFIXER_AVAILABLE = True
except ImportError:
    _PDBFIXER_AVAILABLE = False

# Force obabel — skip Meeko entirely (API changed in newer versions)
_MEEKO_AVAILABLE = False

JOBS_DIR = Path("storage/jobs")


def prepare_receptor(pdb_path: str, job_id: str) -> str:
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    fixed_pdb_path = str(job_dir / "receptor_fixed.pdb")
    pdbqt_path     = str(job_dir / "receptor.pdbqt")

    if _PDBFIXER_AVAILABLE:
        _fix_pdb_with_pdbfixer(pdb_path, fixed_pdb_path)
    else:
        print("[ProteinPrep] PDBFixer not available — using raw PDB.")
        _strip_waters_manually(pdb_path, fixed_pdb_path)

    _convert_pdb_to_pdbqt(fixed_pdb_path, pdbqt_path)
    return pdbqt_path


def prepare_ligand(smiles: str, job_id: str, compound_index: int) -> str:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    job_dir = JOBS_DIR / job_id / "ligands"
    job_dir.mkdir(parents=True, exist_ok=True)

    sdf_path   = str(job_dir / f"ligand_{compound_index}.sdf")
    pdbqt_path = str(job_dir / f"ligand_{compound_index}.pdbqt")

    print(f"[LigandPrep] Preparing compound {compound_index}: {smiles[:40]}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: '{smiles}'")

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        result = AllChem.EmbedMolecule(mol, AllChem.EmbedParameters())
    if result == -1:
        raise RuntimeError(f"RDKit could not generate 3D coords for: '{smiles}'")

    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    print(f"[LigandPrep] SDF written: {sdf_path}")

    # Always use obabel — Meeko API is broken on this env
    _convert_sdf_to_pdbqt_obabel(sdf_path, pdbqt_path)

    print(f"[LigandPrep] PDBQT written: {pdbqt_path}")
    return pdbqt_path


def _fix_pdb_with_pdbfixer(input_pdb: str, output_pdb: str) -> None:
    fixer = PDBFixer(filename=input_pdb)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(output_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


def _strip_waters_manually(input_pdb: str, output_pdb: str) -> None:
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
    cmd = [
        "obabel",
        input_pdb,
        "-O", output_pdbqt,
        "-xr",
        "--partialcharge", "gasteiger",
    ]
    print(f"[ReceptorPrep] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    print(f"[ReceptorPrep] obabel stdout: {result.stdout}")
    print(f"[ReceptorPrep] obabel stderr: {result.stderr}")

    if result.returncode != 0:
        raise RuntimeError(
            f"Open Babel receptor conversion failed:\n"
            f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
        )
    if not Path(output_pdbqt).exists():
        raise RuntimeError(f"obabel ran but no output: {output_pdbqt}")


def _convert_sdf_to_pdbqt_obabel(sdf_path: str, output_pdbqt: str) -> None:
    cmd = [
        "obabel",
        sdf_path,
        "-O", output_pdbqt,
        "--partialcharge", "gasteiger",
        "-xh",
    ]
    print(f"[LigandPrep] Running obabel: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    print(f"[LigandPrep] obabel stdout: {result.stdout}")
    print(f"[LigandPrep] obabel stderr: {result.stderr}")
    print(f"[LigandPrep] obabel returncode: {result.returncode}")

    if result.returncode != 0:
        raise RuntimeError(
            f"Open Babel ligand conversion failed:\n"
            f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
        )
    if not Path(output_pdbqt).exists():
        raise RuntimeError(f"obabel ran but no output: {output_pdbqt}")

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


JOBS_DIR = Path("storage/jobs")


def prepare_receptor(pdb_path: str, job_id: str) -> str:
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    fixed_pdb_path = str(job_dir / "receptor_fixed.pdb")
    pdbqt_path     = str(job_dir / "receptor.pdbqt")

    if _PDBFIXER_AVAILABLE:
        _fix_pdb_with_pdbfixer(pdb_path, fixed_pdb_path)
    else:
        print("[ProteinPrep] PDBFixer not available — stripping waters + HETATM manually.")
        _strip_waters_and_hetatm(pdb_path, fixed_pdb_path)

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


def _strip_waters_and_hetatm(input_pdb: str, output_pdb: str) -> None:
    """
    Strips ALL HETATM records (ligands, cofactors, waters) from PDB.
    Keeps only ATOM records (protein backbone).
    NOTE: Binding site coordinates are computed from raw_pdb_path BEFORE
    this step runs, so stripping here does NOT affect auto-boxing accuracy.
    """
    water_res = {"HOH", "WAT", "DOD"}
    lines_out = []
    with open(input_pdb, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record == "HETATM":
                continue
            if record == "ATOM":
                if line[17:20].strip() in water_res:
                    continue
            lines_out.append(line)

    with open(output_pdb, "w") as f:
        f.writelines(lines_out)

    print(f"[ProteinPrep] Stripped PDB written: {output_pdb} ({len(lines_out)} lines)")


def _convert_pdb_to_pdbqt(input_pdb: str, output_pdbqt: str) -> None:
    cmd = [
        "obabel",
        input_pdb,
        "-e",
        "-O", output_pdbqt,
        "-xr",
        "--partialcharge", "gasteiger",
    ]
    print(f"[ReceptorPrep] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    print(f"[ReceptorPrep] obabel stdout: {result.stdout}")
    print(f"[ReceptorPrep] obabel stderr: {result.stderr}")

    if not Path(output_pdbqt).exists():
        raise RuntimeError(f"obabel produced no receptor file: {output_pdbqt}")

    file_size = Path(output_pdbqt).stat().st_size
    print(f"[ReceptorPrep] Output file size: {file_size} bytes")

    if file_size < 500:
        raise RuntimeError(
            f"Receptor PDBQT too small ({file_size} bytes) — obabel likely failed.\n"
            f"stderr: {result.stderr}"
        )


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
        raise RuntimeError(f"obabel ligand failed:\n{result.stderr}")
    if not Path(output_pdbqt).exists() or Path(output_pdbqt).stat().st_size < 50:
        raise RuntimeError(f"obabel produced empty ligand file: {output_pdbqt}")


def cleanup_docking_files(job_id: str) -> None:
    import shutil
    ligand_dir = JOBS_DIR / job_id / "ligands"
    if ligand_dir.exists():
        shutil.rmtree(ligand_dir)

    out_dir = JOBS_DIR / job_id / "docking_outputs"
    if out_dir.exists():
        for f in out_dir.glob("out_*.pdbqt"):
            f.unlink(missing_ok=True)

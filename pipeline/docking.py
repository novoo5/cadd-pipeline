import subprocess
import os
import re
import asyncio
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from backend.models.schemas import (
    DockingResult,
    DockingPose,
    DockingSpeed,
    BindingSiteMode,
    BindingSiteCoords,
    BindingSiteResidues,
)
from backend.pipeline.protein_prep import prepare_ligand
from backend.utils.binding_site import get_binding_site_box

JOBS_DIR = Path("storage/jobs")

# GNINA exhaustiveness values per speed setting
EXHAUSTIVENESS_MAP = {
    DockingSpeed.fast:     8,
    DockingSpeed.balanced: 16,
    DockingSpeed.thorough: 32,
}

# GNINA docking hard limits
NUM_MODES        = 9    # number of poses to generate per compound
ENERGY_RANGE     = 3.0  # kcal/mol range from best pose to include


def dock_compounds(
    smiles_list: List[str],
    receptor_pdbqt: str,
    job_id: str,
    raw_pdb_path: str,
    binding_site_mode: BindingSiteMode,
    docking_speed: DockingSpeed,
    binding_site_coords: Optional[BindingSiteCoords] = None,
    binding_site_residues: Optional[BindingSiteResidues] = None,
) -> Tuple[List[str], List[DockingResult]]:
    """
    Runs GNINA docking sequentially for each compound in smiles_list.
    Sequential (not parallel) to avoid OOM on free HF 16GB RAM tier.

    Returns:
      - passed_smiles: compounds with best_affinity_kcal < DOCKING_PASS_THRESHOLD
      - all_results: DockingResult for every compound
    """
    # Compute binding site box once (same for all ligands)
    from backend.models.schemas import BindingSiteMode as BSM
    center, box_size = get_binding_site_box(
        pdb_path=raw_pdb_path,
        mode=binding_site_mode,
        coords=binding_site_coords,
        residues=binding_site_residues,
    )
    cx, cy, cz = center
    exhaustiveness = EXHAUSTIVENESS_MAP[docking_speed]

    passed_smiles: List[str] = []
    all_results: List[DockingResult] = []

    for idx, smi in enumerate(smiles_list):
        print(f"[Docking] Compound {idx + 1}/{len(smiles_list)}: {smi[:40]}...")
        try:
            result = _dock_single(
                smiles=smi,
                receptor_pdbqt=receptor_pdbqt,
                job_id=job_id,
                compound_index=idx,
                cx=cx, cy=cy, cz=cz,
                box_size=box_size,
                exhaustiveness=exhaustiveness,
            )
            all_results.append(result)
            if result.passed:
                passed_smiles.append(smi)

        except Exception as e:
            print(f"[Docking] Failed for compound {idx}: {e}")
            all_results.append(DockingResult(
                passed=False,
                best_affinity_kcal=0.0,
                cnn_score=0.0,
                poses=[],
            ))

    return passed_smiles, all_results


def _dock_single(
    smiles: str,
    receptor_pdbqt: str,
    job_id: str,
    compound_index: int,
    cx: float,
    cy: float,
    cz: float,
    box_size: float,
    exhaustiveness: int,
) -> DockingResult:
    """
    Prepares one ligand and runs GNINA docking against the prepared receptor.
    Parses GNINA's output log to extract poses, affinities, and CNN scores.
    """
    # Prepare ligand → .pdbqt
    ligand_pdbqt = prepare_ligand(smiles, job_id, compound_index)

    # Output pose file
    out_dir = JOBS_DIR / job_id / "docking_outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = str(out_dir / f"out_{compound_index}.pdbqt")
    log_file  = str(out_dir / f"log_{compound_index}.txt")

    # Build GNINA command
    # --cnn_scoring rescore: use CNN to rescore Vina poses (best accuracy/speed tradeoff)
    # --cnn crossdock_default2018: pre-trained CNN weights shipped with GNINA
    cmd = [
        "gnina",
        "--receptor",       receptor_pdbqt,
        "--ligand",         ligand_pdbqt,
        "--out",            out_pdbqt,
        "--log",            log_file,
        "--center_x",       str(round(cx, 3)),
        "--center_y",       str(round(cy, 3)),
        "--center_z",       str(round(cz, 3)),
        "--size_x",         str(round(box_size, 1)),
        "--size_y",         str(round(box_size, 1)),
        "--size_z",         str(round(box_size, 1)),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes",      str(NUM_MODES),
        "--energy_range",   str(ENERGY_RANGE),
        "--cnn_scoring",    "rescore",
        "--cnn",            "crossdock_default2018",
        "--cpu",            "2",     # match free HF tier vCPU count
        "--no_gpu",                  # free HF Space has no GPU
    ]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=900,    # 15 min hard timeout per compound
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"GNINA timed out after 15 minutes for compound {compound_index}. "
            "Try using 'Fast' docking speed."
        )

    if proc.returncode != 0:
        raise RuntimeError(
            f"GNINA exited with code {proc.returncode}:\n"
            f"STDERR: {proc.stderr[:500]}"
        )

    # Parse log file for pose results
    poses = _parse_gnina_log(log_file)
    if not poses:
        raise RuntimeError(f"GNINA produced no poses for compound {compound_index}.")

    best_pose = poses[0]
    best_affinity = best_pose.affinity_kcal
    best_cnn_score = _extract_cnn_score(log_file)

    # Pass threshold: affinity must be at least -6.0 kcal/mol
    # (rough rule-of-thumb: below -6 is considered potentially active)
    passed = best_affinity <= -6.0

    return DockingResult(
        passed=passed,
        best_affinity_kcal=round(best_affinity, 4),
        cnn_score=round(best_cnn_score, 4),
        poses=poses,
    )


def _parse_gnina_log(log_file: str) -> List[DockingPose]:
    """
    Parses GNINA's text log file to extract pose table.

    GNINA log format (relevant section):
    -----+------------+----------+----------
    mode | affinity   | CNN score| dist from best mode
         | (kcal/mol) |          | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1 |     -8.412 |    0.742 |    0.000 |    0.000
       2 |     -7.981 |    0.698 |    1.243 |    2.107
    """
    poses = []
    if not Path(log_file).exists():
        return poses

    with open(log_file, "r") as f:
        content = f.read()

    # Regex to match pose table rows
    pattern = re.compile(
        r"^\s*(\d+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)",
        re.MULTILINE,
    )
    for match in pattern.finditer(content):
        try:
            poses.append(DockingPose(
                rank=int(match.group(1)),
                affinity_kcal=float(match.group(2)),
                rmsd_lb=float(match.group(4)),
                rmsd_ub=float(match.group(5)),
            ))
        except (ValueError, IndexError):
            continue

    return poses


def _extract_cnn_score(log_file: str) -> float:
    """
    Extracts the best CNN pose score from GNINA log.
    Returns 0.0 if not found.
    """
    if not Path(log_file).exists():
        return 0.0
    with open(log_file, "r") as f:
        content = f.read()

    # Look for CNN score in the first pose row
    pattern = re.compile(
        r"^\s*1\s*\|\s*[-\d.]+\s*\|\s*([-\d.]+)",
        re.MULTILINE,
    )
    match = pattern.search(content)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return 0.0
    return 0.0


def cleanup_docking_files(job_id: str) -> None:
    """
    Removes intermediate ligand .pdbqt and docking output files after results
    are saved. Keeps log files for debugging but removes large pose files.
    """
    import shutil
    ligand_dir = JOBS_DIR / job_id / "ligands"
    if ligand_dir.exists():
        shutil.rmtree(ligand_dir)

    out_dir = JOBS_DIR / job_id / "docking_outputs"
    if out_dir.exists():
        for f in out_dir.glob("out_*.pdbqt"):
            f.unlink(missing_ok=True)

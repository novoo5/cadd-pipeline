# backend/pipeline/docking.py
import re
import subprocess
import traceback
from pathlib import Path
from typing import List, Tuple, Optional

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

EXHAUSTIVENESS_MAP = {
    DockingSpeed.fast:     4,
    DockingSpeed.balanced: 8,
    DockingSpeed.thorough: 16,
}

NUM_MODES      = 3
ENERGY_RANGE   = 3.0
BOX_SIZE       = 25.0
PASS_THRESHOLD = -5.0


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

    center, _ = get_binding_site_box(
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
                exhaustiveness=exhaustiveness,
            )
            all_results.append(result)
            if result.passed:
                passed_smiles.append(smi)

        except Exception as e:
            print(f"[Docking] FULL ERROR for compound {idx}:")
            print(traceback.format_exc())
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
    exhaustiveness: int,
) -> DockingResult:

    ligand_pdbqt = prepare_ligand(smiles, job_id, compound_index)
    print(f"[Docking] Ligand prepared: {ligand_pdbqt}")

    out_dir = JOBS_DIR / job_id / "docking_outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = str(out_dir / f"out_{compound_index}.pdbqt")

    cmd = [
        "vina",
        "--receptor",       receptor_pdbqt,
        "--ligand",         ligand_pdbqt,
        "--out",            out_pdbqt,
        "--center_x",       str(round(cx, 3)),
        "--center_y",       str(round(cy, 3)),
        "--center_z",       str(round(cz, 3)),
        "--size_x",         str(BOX_SIZE),
        "--size_y",         str(BOX_SIZE),
        "--size_z",         str(BOX_SIZE),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes",      str(NUM_MODES),
        "--energy_range",   str(ENERGY_RANGE),
        "--cpu",            "2",
    ]

    print(f"[Docking] CMD: {' '.join(cmd)}")

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Vina timed out for compound {compound_index}")

    print(f"[Docking] returncode: {proc.returncode}")
    print(f"[Docking] stdout: {proc.stdout[:3000]}")
    print(f"[Docking] stderr: {proc.stderr[:500]}")

    if proc.returncode != 0:
        raise RuntimeError(
            f"Vina failed with code {proc.returncode}:\n{proc.stderr[:500]}"
        )

    poses = _parse_vina_output(proc.stdout)
    print(f"[Docking] Parsed {len(poses)} poses from stdout")

    if not poses:
        print(f"[Docking] WARNING: No poses parsed for compound {compound_index}. Raw stdout tail:")
        print(proc.stdout[-500:])
        return DockingResult(
            passed=False,
            best_affinity_kcal=0.0,
            cnn_score=0.0,
            poses=[],
        )

    best_affinity = poses[0].affinity_kcal
    print(f"[Docking] Best affinity: {best_affinity} kcal/mol")

    passed = best_affinity <= PASS_THRESHOLD

    return DockingResult(
        passed=passed,
        best_affinity_kcal=round(best_affinity, 4),
        cnn_score=0.0,
        poses=poses,
    )


def _parse_vina_output(stdout: str) -> List[DockingPose]:
    """
    Regex-based parser robust against whitespace variations.
    Matches lines like:
       1       -6.642          0          0
       2       -6.479      23.55      25.65
    """
    pattern = re.compile(
        r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s*$"
    )

    poses = []
    in_table = False

    for line in stdout.splitlines():
        if re.match(r"\s*-{5,}\+", line):
            in_table = True
            continue

        if not in_table:
            continue

        m = pattern.match(line)
        if m:
            try:
                poses.append(DockingPose(
                    rank=int(m.group(1)),
                    affinity_kcal=round(float(m.group(2)), 4),
                    rmsd_lb=round(float(m.group(3)), 3),
                    rmsd_ub=round(float(m.group(4)), 3),
                ))
            except Exception as e:
                print(f"[Docking] Parse warning on line '{line.strip()}': {e}")
                continue
        elif in_table and line.strip() == "":
            if poses:
                break

    return poses


def cleanup_docking_files(job_id: str) -> None:
    import shutil
    ligand_dir = JOBS_DIR / job_id / "ligands"
    if ligand_dir.exists():
        shutil.rmtree(ligand_dir)

    out_dir = JOBS_DIR / job_id / "docking_outputs"
    if out_dir.exists():
        for f in out_dir.glob("out_*.pdbqt"):
            f.unlink(missing_ok=True)

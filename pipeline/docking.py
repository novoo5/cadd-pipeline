# backend/pipeline/docking.py
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

NUM_MODES    = 3
ENERGY_RANGE = 3.0


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


def _get_ligand_box_size(ligand_pdbqt: str, padding: float = 8.0) -> float:
    """
    Calculates a dynamic docking box size from the ligand's actual 3D
    atom coordinates + padding. Clamps between 20 and 30 Å.
    Prevents grid overflow for large molecules and wasted space for small ones.
    """
    coords = []
    try:
        with open(ligand_pdbqt) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
        if not coords:
            return 25.0
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]
        span = max(
            max(xs) - min(xs),
            max(ys) - min(ys),
            max(zs) - min(zs),
        )
        size = max(20.0, min(span + padding, 30.0))
        return round(size, 1)
    except Exception:
        return 25.0


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

    # ✅ Dynamic box size from actual ligand 3D dimensions
    dynamic_box = _get_ligand_box_size(ligand_pdbqt)
    print(f"[Docking] Dynamic box size: {dynamic_box} Å")

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
        "--size_x",         str(dynamic_box),
        "--size_y",         str(dynamic_box),
        "--size_z",         str(dynamic_box),
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
    print(f"[Docking] stdout: {proc.stdout[:2000]}")
    print(f"[Docking] stderr: {proc.stderr[:500]}")

    if proc.returncode != 0:
        raise RuntimeError(
            f"Vina failed with code {proc.returncode}:\n{proc.stderr[:500]}"
        )

    poses = _parse_vina_output(proc.stdout)
    if not poses:
        raise RuntimeError(f"Vina produced no poses for compound {compound_index}")

    best_affinity = poses[0].affinity_kcal
    print(f"[Docking] Best affinity: {best_affinity} kcal/mol")

    passed = best_affinity <= -6.0

    return DockingResult(
        passed=passed,
        best_affinity_kcal=round(best_affinity, 4),
        cnn_score=0.0,
        poses=poses,
    )


def _parse_vina_output(stdout: str) -> List[DockingPose]:
    """
    Parses AutoDock Vina stdout table:
    -----+------------+----------+----------
       1       -8.412      0.000      0.000
    """
    poses = []
    in_table = False
    for line in stdout.splitlines():
        if line.strip().startswith("-----+"):
            in_table = True
            continue
        if in_table:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    mode     = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb  = float(parts[2])
                    rmsd_ub  = float(parts[3])
                    poses.append(DockingPose(
                        pose_number=mode,
                        affinity_kcal=round(affinity, 4),
                        rmsd_lower=round(rmsd_lb, 3),
                        rmsd_upper=round(rmsd_ub, 3),
                    ))
                except (ValueError, IndexError):
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

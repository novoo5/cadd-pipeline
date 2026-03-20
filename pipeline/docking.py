# backend/pipeline/docking.py
import os
import re
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
                box_size=box_size,
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
    box_size: float,
    exhaustiveness: int,
) -> DockingResult:
    from vina import Vina

    ligand_pdbqt = prepare_ligand(smiles, job_id, compound_index)
    print(f"[Docking] Ligand prepared: {ligand_pdbqt}")

    out_dir = JOBS_DIR / job_id / "docking_outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = str(out_dir / f"out_{compound_index}.pdbqt")

    # Initialize AutoDock Vina
    v = Vina(sf_name='vina', cpu=2, verbosity=1)

    print(f"[Docking] Loading receptor: {receptor_pdbqt}")
    v.set_receptor(receptor_pdbqt)

    print(f"[Docking] Loading ligand: {ligand_pdbqt}")
    v.set_ligand_from_file(ligand_pdbqt)

    # Define search box
    v.compute_vina_maps(
        center=[cx, cy, cz],
        box_size=[box_size, box_size, box_size],
    )

    print(f"[Docking] Running docking with exhaustiveness={exhaustiveness}...")
    v.dock(exhaustiveness=exhaustiveness, n_poses=NUM_MODES)

    # Get scores
    energies = v.energies(n_poses=NUM_MODES)
    print(f"[Docking] Raw energies: {energies}")

    if energies is None or len(energies) == 0:
        raise RuntimeError(f"Vina returned no poses for compound {compound_index}")

    # Write output
    v.write_poses(out_pdbqt, n_poses=NUM_MODES, overwrite=True)

    poses = []
    for i, energy in enumerate(energies):
        try:
            poses.append(DockingPose(
                rank=i + 1,
                affinity_kcal=float(energy[0]),  # index 0 = total energy
                rmsd_lb=float(energy[1]) if len(energy) > 1 else 0.0,
                rmsd_ub=float(energy[2]) if len(energy) > 2 else 0.0,
            ))
        except (ValueError, IndexError) as e:
            print(f"[Docking] Pose parse error at index {i}: {e}")
            continue

    if not poses:
        raise RuntimeError(f"Could not parse any poses for compound {compound_index}")

    best_affinity = poses[0].affinity_kcal
    print(f"[Docking] Best affinity: {best_affinity} kcal/mol")

    passed = best_affinity <= -6.0

    return DockingResult(
        passed=passed,
        best_affinity_kcal=round(best_affinity, 4),
        cnn_score=0.0,  # Vina doesn't have CNN score, that's GNINA-specific
        poses=poses,
    )


def cleanup_docking_files(job_id: str) -> None:
    import shutil
    ligand_dir = JOBS_DIR / job_id / "ligands"
    if ligand_dir.exists():
        shutil.rmtree(ligand_dir)

    out_dir = JOBS_DIR / job_id / "docking_outputs"
    if out_dir.exists():
        for f in out_dir.glob("out_*.pdbqt"):
            f.unlink(missing_ok=True)

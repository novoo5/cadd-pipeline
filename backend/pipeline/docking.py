"""
Molecular docking via AutoDock Vina.

New in this version (download feature):
  - _dock_single() now accepts raw_pdb_path and generates two extra files after docking:
      best_pose_{idx}_H.pdb           — Pose 1 only, with hydrogens added (obabel)
      complex_{idx}_PLIP_ready.pdb    — receptor ATOM records + docked ligand HETATM,
                                        ready for PLIP / PyMOL without any CLI work
  - cleanup_docking_files() now only removes the ligand preparation directory.
    docking_outputs/ is intentionally preserved so users can download all three files.
"""

import re
import subprocess
import traceback
from pathlib import Path
from typing import List, Tuple, Optional, Callable

from backend.models.schemas import (
    DockingResult,
    DockingPose,
    DockingSpeed,
    BindingSiteMode,
    BindingSiteCoords,
    BindingSiteResidues,
    BindingSiteInfo,
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


# ── Public API ────────────────────────────────────────────────────────────────

def dock_compounds(
    smiles_list:           List[str],
    receptor_pdbqt:        str,
    job_id:                str,
    raw_pdb_path:          str,
    binding_site_mode:     BindingSiteMode,
    docking_speed:         DockingSpeed,
    binding_site_coords:   Optional[BindingSiteCoords]   = None,
    binding_site_residues: Optional[BindingSiteResidues] = None,
    on_progress:           Optional[Callable[[int, int], None]] = None,
) -> Tuple[List[str], List[DockingResult], BindingSiteInfo]:

    center, _, binding_site_info = get_binding_site_box(
        pdb_path=raw_pdb_path,
        mode=binding_site_mode,
        coords=binding_site_coords,
        residues=binding_site_residues,
    )
    cx, cy, cz = center
    exhaustiveness = EXHAUSTIVENESS_MAP[docking_speed]

    passed_smiles: List[str]        = []
    all_results:   List[DockingResult] = []
    total = len(smiles_list)

    for idx, smi in enumerate(smiles_list):
        print(f"[Docking] Compound {idx + 1}/{total}: {smi[:40]}...")
        try:
            result = _dock_single(
                smiles=smi,
                receptor_pdbqt=receptor_pdbqt,
                raw_pdb_path=raw_pdb_path,   # ← passed through for complex generation
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

        if on_progress:
            on_progress(idx + 1, total)

    return passed_smiles, all_results, binding_site_info


# ── Single-compound docking ───────────────────────────────────────────────────

def _dock_single(
    smiles:          str,
    receptor_pdbqt:  str,
    raw_pdb_path:    str,   # ← NEW: needed for complex PDB generation
    job_id:          str,
    compound_index:  int,
    cx:              float,
    cy:              float,
    cz:              float,
    exhaustiveness:  int,
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
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"Vina timed out for compound {compound_index}")

    print(f"[Docking] returncode: {proc.returncode}")
    print(f"[Docking] stdout: {proc.stdout[:3000]}")
    print(f"[Docking] stderr: {proc.stderr[:500]}")

    if proc.returncode != 0:
        raise RuntimeError(f"Vina failed with code {proc.returncode}:\n{proc.stderr[:500]}")

    poses = _parse_vina_output(proc.stdout)
    print(f"[Docking] Parsed {len(poses)} poses from stdout")

    if not poses:
        print(f"[Docking] WARNING: No poses for compound {compound_index}.")
        print(proc.stdout[-500:])
        return DockingResult(passed=False, best_affinity_kcal=0.0, cnn_score=0.0, poses=[])

    best_affinity = poses[0].affinity_kcal
    print(f"[Docking] Best affinity: {best_affinity} kcal/mol")

    # ── Post-docking: generate downloadable PDB files ────────────────────────
    pose_pdb_path = out_dir / f"best_pose_{compound_index}_H.pdb"
    _generate_pose_pdb(out_pdbqt, pose_pdb_path)
    _generate_complex_pdb(
        pose_pdb_path=pose_pdb_path,
        receptor_pdb_path=raw_pdb_path,
        out_path=out_dir / f"complex_{compound_index}_PLIP_ready.pdb",
        compound_index=compound_index,
    )

    return DockingResult(
        passed=best_affinity <= PASS_THRESHOLD,
        best_affinity_kcal=round(best_affinity, 4),
        cnn_score=0.0,
        poses=poses,
    )


# ── Post-docking file generators ──────────────────────────────────────────────

def _generate_pose_pdb(out_pdbqt: str, pose_pdb_path: Path) -> None:
    """
    Use obabel to extract pose 1 from out_{idx}.pdbqt and add hydrogens.
    Output: best_pose_{idx}_H.pdb — ready for PLIP or PyMOL.
    Silently skips if obabel is not installed or fails.
    """
    cmd = [
        "obabel", out_pdbqt,
        "-O", str(pose_pdb_path),
        "-f", "1", "-l", "1",   # first MODEL block only
        "-h",                    # add hydrogens
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            print(f"[Docking] obabel pose→pdb warning: {result.stderr[:300]}")
        else:
            print(f"[Docking] Pose PDB written: {pose_pdb_path}")
    except FileNotFoundError:
        print("[Docking] obabel not found — skipping pose PDB generation.")
    except subprocess.TimeoutExpired:
        print(f"[Docking] obabel timed out for {pose_pdb_path.name}.")


def _generate_complex_pdb(
    pose_pdb_path:    Path,
    receptor_pdb_path: str,
    out_path:         Path,
    compound_index:   int,
) -> None:
    """
    Merge receptor + docked ligand into a single PDB file for PLIP.

    Receptor handling:
      - ATOM lines are always kept (protein backbone + sidechains).
      - HETATM lines are kept only for single-element entries (structural ions
        like Zn²⁺, Mg²⁺, Ca²⁺) — crystal ligands are stripped to avoid PLIP
        confusion.
      - END / CONECT / REMARK from the receptor are dropped.

    Ligand handling:
      - All ATOM/HETATM lines from the pose PDB are re-labelled as HETATM so
        PLIP correctly identifies them as the bound ligand.
      - CONECT records are preserved for bond topology.

    Skips silently if pose_pdb_path does not exist (obabel unavailable).
    """
    if not pose_pdb_path.exists():
        print(f"[Docking] Pose PDB missing — skipping complex generation for compound {compound_index}.")
        return

    try:
        receptor_lines = Path(receptor_pdb_path).read_text(errors="ignore").splitlines()
        pose_lines     = pose_pdb_path.read_text(errors="ignore").splitlines()

        receptor_block: List[str] = []
        for line in receptor_lines:
            rec = line[:6].strip()
            if rec == "ATOM":
                receptor_block.append(line)
            elif rec == "HETATM":
                # Keep only single-atom HETATM (structural metal ions)
                atom_name = line[12:16].strip()
                if len(atom_name) <= 2:
                    receptor_block.append(line)

        ligand_block: List[str] = []
        for line in pose_lines:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                # Force HETATM label so PLIP treats this as the ligand
                ligand_block.append("HETATM" + line[6:])
            elif rec == "CONECT":
                ligand_block.append(line)

        combined_lines = receptor_block + ["TER"] + ligand_block + ["END"]
        out_path.write_text("\n".join(combined_lines) + "\n")
        print(f"[Docking] Complex PDB written: {out_path}")

    except Exception as e:
        print(f"[Docking] Complex generation error for compound {compound_index}: {e}")


# ── Vina stdout parser ────────────────────────────────────────────────────────

def _parse_vina_output(stdout: str) -> List[DockingPose]:
    pattern = re.compile(
        r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s*$"
    )
    poses    = []
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


# ── Cleanup ───────────────────────────────────────────────────────────────────

def cleanup_docking_files(job_id: str) -> None:
    """
    Removes ligand *preparation* files only (input .pdbqt files generated
    by prepare_ligand).

    docking_outputs/ is intentionally preserved — it contains the three
    downloadable files:
        out_{idx}.pdbqt                 — raw Vina output
        best_pose_{idx}_H.pdb           — best pose with hydrogens
        complex_{idx}_PLIP_ready.pdb    — receptor + ligand complex
    """
    import shutil
    ligand_dir = JOBS_DIR / job_id / "ligands"
    if ligand_dir.exists():
        shutil.rmtree(ligand_dir)
    # NOTE: docking_outputs/ is NOT cleaned — files are needed for user download.

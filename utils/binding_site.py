import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List
from backend.models.schemas import (
    BindingSiteMode,
    BindingSiteCoords,
    BindingSiteResidues,
)


def get_binding_site_box(
    pdb_path: str,
    mode: BindingSiteMode,
    coords: Optional[BindingSiteCoords] = None,
    residues: Optional[BindingSiteResidues] = None,
) -> Tuple[Tuple[float, float, float], float]:
    """
    Master function — returns ((cx, cy, cz), box_size) for GNINA docking.

    Dispatches to the correct method based on mode:
      - auto        → detects binding site from HETATM ligand centroid, or uses full protein centroid
      - coordinates → uses user-supplied centroid + box size directly
      - residues    → computes centroid from specified residue range
    """
    if mode == BindingSiteMode.coordinates:
        assert coords is not None
        return (coords.x, coords.y, coords.z), coords.box_size

    if mode == BindingSiteMode.residues:
        assert residues is not None
        return _centroid_from_residues(pdb_path, residues.chain, residues.residue_start, residues.residue_end)

    # Auto mode
    return _auto_detect_binding_site(pdb_path)


def _parse_atom_coords(pdb_path: str, record_types: List[str], chain: Optional[str] = None) -> np.ndarray:
    """
    Parses PDB ATOM/HETATM records and returns an (N, 3) numpy array of coordinates.
    Optionally filters by chain ID.
    Skips water molecules (HOH, WAT).
    """
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record not in record_types:
                continue

            # Skip water
            res_name = line[17:20].strip()
            if res_name in ("HOH", "WAT", "DOD"):
                continue

            if chain is not None:
                chain_id = line[21].strip()
                if chain_id != chain:
                    continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except ValueError:
                continue

    if not coords:
        return np.array([]).reshape(0, 3)
    return np.array(coords)


def _auto_detect_binding_site(pdb_path: str) -> Tuple[Tuple[float, float, float], float]:
    """
    Strategy 1: If a co-crystallized ligand (HETATM) exists, use its centroid as binding site.
    Strategy 2: Fall back to full protein centroid (blind docking, larger box).
    """
    hetatm_coords = _parse_atom_coords(pdb_path, ["HETATM"])

    if len(hetatm_coords) >= 3:
        centroid = hetatm_coords.mean(axis=0)
        # Box size = ligand extent + 10Å padding on each side
        extents = hetatm_coords.max(axis=0) - hetatm_coords.min(axis=0)
        box_size = float(np.max(extents)) + 20.0
        box_size = min(max(box_size, 20.0), 40.0)  # clamp between 20–40Å
        return tuple(centroid.tolist()), box_size  # type: ignore

    # Fallback: full protein centroid, large blind docking box
    atom_coords = _parse_atom_coords(pdb_path, ["ATOM"])
    if len(atom_coords) == 0:
        raise ValueError(
            "No ATOM or HETATM records found in the PDB file. "
            "The file may be corrupted or empty."
        )

    centroid = atom_coords.mean(axis=0)
    extents = atom_coords.max(axis=0) - atom_coords.min(axis=0)
    box_size = float(np.max(extents)) + 10.0
    box_size = min(box_size, 60.0)  # cap at 60Å for sanity

    return tuple(centroid.tolist()), box_size  # type: ignore


def _centroid_from_residues(
    pdb_path: str,
    chain: str,
    residue_start: int,
    residue_end: int,
) -> Tuple[Tuple[float, float, float], float]:
    """
    Computes centroid from all atoms in the given residue range on the given chain.
    Returns ((cx, cy, cz), box_size) with box_size = extent + 12Å padding.
    """
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            chain_id = line[21].strip()
            if chain_id != chain:
                continue
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                continue
            if not (residue_start <= res_seq <= residue_end):
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append([x, y, z])
            except ValueError:
                continue

    if not coords:
        raise ValueError(
            f"No atoms found in chain '{chain}', residues {residue_start}–{residue_end}. "
            "Check chain ID and residue numbers against the PDB file."
        )

    arr = np.array(coords)
    centroid = arr.mean(axis=0)
    extents = arr.max(axis=0) - arr.min(axis=0)
    box_size = float(np.max(extents)) + 12.0
    box_size = min(max(box_size, 15.0), 45.0)

    return tuple(centroid.tolist()), box_size  # type: ignore

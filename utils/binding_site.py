import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List, Dict
from backend.models.schemas import (
    BindingSiteMode,
    BindingSiteCoords,
    BindingSiteResidues,
)


# ── Residues to exclude from native ligand detection ────────────────────────
# These are crystallographic artifacts, buffer molecules, and ions that
# appear as HETATM in PDB files but are NOT the native drug/ligand.
EXCLUDED_HETATM = {
    # Water
    "HOH", "WAT", "DOD", "H2O",
    # Glycols and cryoprotectants
    "EDO", "EOH", "GOL", "PEG", "PGE", "PG4", "MPD", "1PE", "P6G",
    # Organic solvents
    "DMS", "ACT", "ACE", "FMT", "IPA", "MET", "ETA",
    # Buffer molecules
    "TRS", "MES", "HEP", "BIS", "IMD", "EPE", "BME", "DTT", "TCE",
    # Sulfate/phosphate/citrate
    "SO4", "PO4", "SUL", "CIT", "TLA", "FLC", "NO3", "ACT",
    # Metal ions (single atom HETATM)
    "ZN", "MG", "CA", "NA", "K", "MN", "FE", "CU", "NI", "CO",
    "CD", "HG", "PT", "AU", "AG", "IOD", "BR", "CL", "ION",
    # Misc crystallography
    "UNK", "UNL", "LIG", "INH",
}


def get_binding_site_box(
    pdb_path: str,
    mode: BindingSiteMode,
    coords: Optional[BindingSiteCoords] = None,
    residues: Optional[BindingSiteResidues] = None,
) -> Tuple[Tuple[float, float, float], float]:
    """
    Returns ((cx, cy, cz), box_size) for AutoDock Vina docking.

    Modes:
      - auto        → detects native ligand centroid from HETATM records
      - coordinates → uses user-supplied centroid directly
      - residues    → computes centroid from specified residue range
    """
    if mode == BindingSiteMode.coordinates:
        assert coords is not None
        return (coords.x, coords.y, coords.z), coords.box_size

    if mode == BindingSiteMode.residues:
        assert residues is not None
        return _centroid_from_residues(
            pdb_path, residues.chain, residues.residue_start, residues.residue_end
        )

    return _auto_detect_binding_site(pdb_path)


def _auto_detect_binding_site(pdb_path: str) -> Tuple[Tuple[float, float, float], float]:
    """
    Auto-boxing algorithm:
    1. Parse all HETATM records, grouped by (residue_name, chain, residue_seq).
    2. Exclude known buffers/ions/solvents.
    3. Select the group with the MOST atoms = most likely the native co-crystallized ligand.
    4. Calculate centroid of that group as the docking box center.
    5. Fall back to full protein centroid if no valid ligand found.
    """
    groups = _group_hetatm_by_residue(pdb_path)

    # Filter out excluded molecules
    drug_candidates = {
        k: v for k, v in groups.items()
        if k[0] not in EXCLUDED_HETATM and len(v) >= 5  # at least 5 atoms = real molecule
    }

    if drug_candidates:
        # Pick the residue group with the most atoms (= native ligand)
        best_key = max(drug_candidates, key=lambda k: len(drug_candidates[k]))
        native_coords = np.array(drug_candidates[best_key])

        print(f"[BindingSite] Native ligand detected: '{best_key[0]}' "
              f"chain={best_key[1]} res={best_key[2]} "
              f"({len(native_coords)} atoms)")

        centroid = native_coords.mean(axis=0)
        extents = native_coords.max(axis=0) - native_coords.min(axis=0)
        box_size = float(np.max(extents)) + 20.0
        box_size = min(max(box_size, 22.0), 40.0)

        print(f"[BindingSite] Center: X={centroid[0]:.3f} Y={centroid[1]:.3f} Z={centroid[2]:.3f} | Box={box_size:.1f}Å")
        return tuple(centroid.tolist()), box_size  # type: ignore

    # Fallback: if no native ligand, use full protein centroid (blind docking)
    print("[BindingSite] No native ligand found — falling back to protein centroid (blind docking).")
    atom_coords = _parse_atom_coords(pdb_path, ["ATOM"])

    if len(atom_coords) == 0:
        raise ValueError(
            "No ATOM or HETATM records found in the PDB file. "
            "The file may be corrupted or empty."
        )

    centroid = atom_coords.mean(axis=0)
    extents = atom_coords.max(axis=0) - atom_coords.min(axis=0)
    box_size = float(np.max(extents)) + 10.0
    box_size = min(box_size, 60.0)

    print(f"[BindingSite] Fallback center: X={centroid[0]:.3f} Y={centroid[1]:.3f} Z={centroid[2]:.3f} | Box={box_size:.1f}Å")
    return tuple(centroid.tolist()), box_size  # type: ignore


def _group_hetatm_by_residue(pdb_path: str) -> Dict[Tuple[str, str, int], List[List[float]]]:
    """
    Reads all HETATM records and groups atom coordinates by
    (residue_name, chain_id, residue_seq_num).
    Returns a dict of {(res_name, chain, seq): [[x,y,z], ...]}
    """
    groups: Dict[Tuple[str, str, int], List[List[float]]] = {}

    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            try:
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_seq  = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue

            key = (res_name, chain_id, res_seq)
            if key not in groups:
                groups[key] = []
            groups[key].append([x, y, z])

    return groups


def _parse_atom_coords(pdb_path: str, record_types: List[str], chain: Optional[str] = None) -> np.ndarray:
    """
    Parses PDB ATOM/HETATM records and returns an (N, 3) numpy array.
    Skips water molecules.
    """
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record not in record_types:
                continue

            res_name = line[17:20].strip()
            if res_name in ("HOH", "WAT", "DOD"):
                continue

            if chain is not None:
                if line[21].strip() != chain:
                    continue

            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue

    return np.array(coords) if coords else np.array([]).reshape(0, 3)


def _centroid_from_residues(
    pdb_path: str,
    chain: str,
    residue_start: int,
    residue_end: int,
) -> Tuple[Tuple[float, float, float], float]:
    """
    Computes centroid from all atoms in a given residue range on a given chain.
    """
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            if line[21].strip() != chain:
                continue
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                continue
            if not (residue_start <= res_seq <= residue_end):
                continue
            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
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

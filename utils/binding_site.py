import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List, Dict
from backend.models.schemas import (
    BindingSiteMode,
    BindingSiteCoords,
    BindingSiteResidues,
    BindingSiteInfo,
)

EXCLUDED_HETATM = {
    "HOH", "WAT", "DOD", "H2O",
    "EDO", "EOH", "GOL", "PEG", "PGE", "PG4", "MPD", "1PE", "P6G",
    "DMS", "ACT", "ACE", "FMT", "IPA", "MET", "ETA",
    "TRS", "MES", "HEP", "BIS", "IMD", "EPE", "BME", "DTT", "TCE",
    "SO4", "PO4", "SUL", "CIT", "TLA", "FLC", "NO3",
    "ZN", "MG", "CA", "NA", "K", "MN", "FE", "CU", "NI", "CO",
    "CD", "HG", "PT", "AU", "AG", "IOD", "BR", "CL", "ION",
    "UNK", "UNL", "LIG", "INH",
}


def get_binding_site_box(
    pdb_path: str,
    mode: BindingSiteMode,
    coords: Optional[BindingSiteCoords] = None,
    residues: Optional[BindingSiteResidues] = None,
) -> Tuple[Tuple[float, float, float], float, BindingSiteInfo]:
    """
    Returns ((cx, cy, cz), box_size, BindingSiteInfo).
    Third value carries human-readable info about what was detected.
    """
    if mode == BindingSiteMode.coordinates:
        assert coords is not None
        info = BindingSiteInfo(
            detection_mode="user_coordinates",
            detected_ligand_name=None,
            center_x=coords.x, center_y=coords.y, center_z=coords.z,
            box_size=coords.box_size,
        )
        return (coords.x, coords.y, coords.z), coords.box_size, info

    if mode == BindingSiteMode.residues:
        assert residues is not None
        center, box_size = _centroid_from_residues(
            pdb_path, residues.chain, residues.residue_start, residues.residue_end
        )
        cx, cy, cz = center
        info = BindingSiteInfo(
            detection_mode="user_residues",
            detected_ligand_name=f"Chain {residues.chain} res {residues.residue_start}–{residues.residue_end}",
            center_x=round(cx, 3), center_y=round(cy, 3), center_z=round(cz, 3),
            box_size=round(box_size, 1),
        )
        return center, box_size, info

    return _auto_detect_binding_site(pdb_path)


def _auto_detect_binding_site(
    pdb_path: str,
) -> Tuple[Tuple[float, float, float], float, BindingSiteInfo]:
    groups = _group_hetatm_by_residue(pdb_path)

    drug_candidates = {
        k: v for k, v in groups.items()
        if k[0] not in EXCLUDED_HETATM and len(v) >= 5
    }

    if drug_candidates:
        best_key = max(drug_candidates, key=lambda k: len(drug_candidates[k]))
        native_coords = np.array(drug_candidates[best_key])
        ligand_name = best_key[0]

        print(f"[BindingSite] Native ligand detected: '{ligand_name}' "
              f"chain={best_key[1]} res={best_key[2]} ({len(native_coords)} atoms)")

        centroid = native_coords.mean(axis=0)
        extents = native_coords.max(axis=0) - native_coords.min(axis=0)
        box_size = float(np.max(extents)) + 20.0
        box_size = round(min(max(box_size, 22.0), 40.0), 1)
        cx, cy, cz = (round(float(centroid[0]), 3),
                      round(float(centroid[1]), 3),
                      round(float(centroid[2]), 3))

        print(f"[BindingSite] Center: X={cx} Y={cy} Z={cz} | Box={box_size}Å")

        info = BindingSiteInfo(
            detection_mode="native_ligand",
            detected_ligand_name=ligand_name,
            center_x=cx, center_y=cy, center_z=cz,
            box_size=box_size,
        )
        return (cx, cy, cz), box_size, info

    # Fallback: full protein centroid
    print("[BindingSite] No native ligand found — falling back to protein centroid.")
    atom_coords = _parse_atom_coords(pdb_path, ["ATOM"])

    if len(atom_coords) == 0:
        raise ValueError("No ATOM or HETATM records found. PDB file may be corrupted.")

    centroid = atom_coords.mean(axis=0)
    extents = atom_coords.max(axis=0) - atom_coords.min(axis=0)
    box_size = round(min(float(np.max(extents)) + 10.0, 60.0), 1)
    cx, cy, cz = (round(float(centroid[0]), 3),
                  round(float(centroid[1]), 3),
                  round(float(centroid[2]), 3))

    print(f"[BindingSite] Fallback center: X={cx} Y={cy} Z={cz} | Box={box_size}Å")

    info = BindingSiteInfo(
        detection_mode="protein_centroid",
        detected_ligand_name=None,
        center_x=cx, center_y=cy, center_z=cz,
        box_size=box_size,
    )
    return (cx, cy, cz), box_size, info


def _group_hetatm_by_residue(pdb_path: str) -> Dict[Tuple[str, str, int], List[List[float]]]:
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
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            record = line[:6].strip()
            if record not in record_types:
                continue
            if line[17:20].strip() in ("HOH", "WAT", "DOD"):
                continue
            if chain is not None and line[21].strip() != chain:
                continue
            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
    return np.array(coords) if coords else np.array([]).reshape(0, 3)


def _centroid_from_residues(
    pdb_path: str, chain: str, residue_start: int, residue_end: int,
) -> Tuple[Tuple[float, float, float], float]:
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            if line[21].strip() != chain:
                continue
            try:
                if not (residue_start <= int(line[22:26].strip()) <= residue_end):
                    continue
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
    if not coords:
        raise ValueError(f"No atoms in chain '{chain}' res {residue_start}–{residue_end}.")
    arr = np.array(coords)
    centroid = arr.mean(axis=0)
    box_size = round(min(max(float(np.max(arr.max(0) - arr.min(0))) + 12.0, 15.0), 45.0), 1)
    return tuple(centroid.tolist()), box_size  # type: ignore

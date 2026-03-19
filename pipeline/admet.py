import asyncio
from typing import List, Tuple, Dict, Any
from backend.models.schemas import ADMETResult

# ── ADMET-AI import with graceful fallback ───────────────────────────────────
# ADMET-AI uses a Chemprop-RDKit ensemble trained on 41 TDC ADMET datasets.
# pip install admet-ai
try:
    from admet_ai import ADMETModel
    _ADMET_MODEL = ADMETModel()
    _ADMET_AVAILABLE = True
except ImportError:
    _ADMET_AVAILABLE = False
    _ADMET_MODEL = None


# ── Hard-fail thresholds ─────────────────────────────────────────────────────
# These values are conservative clinical cutoffs. Any compound breaching
# these is flagged and fails the ADMET filter.
THRESHOLDS: Dict[str, Any] = {
    "herg_inhibition":       {"max": 0.5,  "label": "High hERG cardiotoxicity risk"},
    "hepatotoxicity":        {"max": 0.5,  "label": "Hepatotoxic"},
    "caco2_permeability":    {"min": -5.15,"label": "Low intestinal permeability"},   # log cm/s
    "oral_bioavailability":  {"min": 0.3,  "label": "Poor oral bioavailability (<30%)"},
}

# BBB penetration is informational only (not a hard fail) — depends on target
# Some CNS drugs need BBB penetration; most non-CNS drugs should not cross it.


def run_admet_batch(smiles_list: List[str]) -> Tuple[List[str], List[ADMETResult]]:
    """
    Runs ADMET-AI prediction on all compounds in batch.
    ADMET-AI processes the full list at once (vectorized inference) — fast.

    Returns:
      - passed_smiles: compounds with no hard-fail flags
      - all_results: ADMETResult for every input compound
    """
    if not _ADMET_AVAILABLE:
        return _fallback_admet(smiles_list)

    try:
        # ADMET-AI accepts a list of SMILES and returns a DataFrame
        df = _ADMET_MODEL.predict(smiles=smiles_list)
        return _parse_admet_dataframe(df, smiles_list)

    except Exception as e:
        print(f"[ADMET] ADMET-AI inference failed: {e}. Using fallback.")
        return _fallback_admet(smiles_list)


def _parse_admet_dataframe(df, smiles_list: List[str]) -> Tuple[List[str], List[ADMETResult]]:
    """
    Parses the ADMET-AI output DataFrame into ADMETResult objects.
    Applies hard-fail thresholds and collects flag strings.
    """
    passed_smiles: List[str] = []
    all_results: List[ADMETResult] = []

    # ADMET-AI column name mappings (TDC dataset names → our schema fields)
    # These are the actual column names returned by ADMET-AI v1.x
    COL_MAP = {
        "herg_inhibition":      "hERG",
        "caco2_permeability":   "Caco2_Wang",
        "bbb_penetration":      "BBB_Martins",
        "hepatotoxicity":       "DILI",
        "oral_bioavailability": "Bioavailability_Ma",
    }

    for i, smi in enumerate(smiles_list):
        try:
            row = df.iloc[i]

            herg   = _safe_get(row, COL_MAP["herg_inhibition"],      0.5)
            caco2  = _safe_get(row, COL_MAP["caco2_permeability"],   -5.5)
            bbb    = _safe_get(row, COL_MAP["bbb_penetration"],       0.5)
            hepato = _safe_get(row, COL_MAP["hepatotoxicity"],        0.5)
            bioav  = _safe_get(row, COL_MAP["oral_bioavailability"],  0.3)

            flags = _compute_flags(herg, caco2, hepato, bioav)
            passed = len(flags) == 0

            result = ADMETResult(
                passed=passed,
                herg_inhibition=round(herg, 4),
                caco2_permeability=round(caco2, 4),
                bbb_penetration=round(bbb, 4),
                hepatotoxicity=round(hepato, 4),
                oral_bioavailability=round(bioav, 4),
                flags=flags,
            )
            all_results.append(result)
            if passed:
                passed_smiles.append(smi)

        except Exception as e:
            print(f"[ADMET] Error parsing row {i} for SMILES '{smi}': {e}")
            all_results.append(_make_failed_result(["ADMET parse error"]))

    return passed_smiles, all_results


def _compute_flags(herg: float, caco2: float, hepato: float, bioav: float) -> List[str]:
    """Applies hard-fail threshold rules and returns a list of flag strings."""
    flags = []
    if herg   > THRESHOLDS["herg_inhibition"]["max"]:
        flags.append(THRESHOLDS["herg_inhibition"]["label"])
    if hepato > THRESHOLDS["hepatotoxicity"]["max"]:
        flags.append(THRESHOLDS["hepatotoxicity"]["label"])
    if caco2  < THRESHOLDS["caco2_permeability"]["min"]:
        flags.append(THRESHOLDS["caco2_permeability"]["label"])
    if bioav  < THRESHOLDS["oral_bioavailability"]["min"]:
        flags.append(THRESHOLDS["oral_bioavailability"]["label"])
    return flags


def _safe_get(row, col_name: str, default: float) -> float:
    """Safely extracts a float from a DataFrame row, with a default fallback."""
    try:
        val = row.get(col_name, default)
        return float(val) if val is not None else default
    except (TypeError, ValueError):
        return default


def _make_failed_result(flags: List[str]) -> ADMETResult:
    return ADMETResult(
        passed=False,
        herg_inhibition=1.0,
        caco2_permeability=-6.0,
        bbb_penetration=0.0,
        hepatotoxicity=1.0,
        oral_bioavailability=0.0,
        flags=flags,
    )


def _fallback_admet(smiles_list: List[str]) -> Tuple[List[str], List[ADMETResult]]:
    """
    Fallback when ADMET-AI is not installed.
    Uses a simple RDKit-based rule-of-thumb heuristic (much less accurate).
    All compounds pass with a warning flag so the pipeline can still run.
    This should NOT be used in production — install admet-ai properly.
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    results = []
    passed = []

    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                results.append(_make_failed_result(["Invalid SMILES in fallback"]))
                continue

            # Rough hERG estimate: high logP + aromatic rings = higher risk
            logp = Descriptors.MolLogP(mol)
            herg_rough = min(0.1 + max(0.0, (logp - 2.0) * 0.1), 1.0)

            result = ADMETResult(
                passed=True,  # permissive in fallback
                herg_inhibition=round(herg_rough, 4),
                caco2_permeability=-4.5,
                bbb_penetration=0.3,
                hepatotoxicity=0.3,
                oral_bioavailability=0.5,
                flags=["⚠️ ADMET-AI not installed — using rough heuristic estimates"],
            )
            results.append(result)
            passed.append(smi)

        except Exception:
            results.append(_make_failed_result(["Fallback ADMET error"]))

    return passed, results

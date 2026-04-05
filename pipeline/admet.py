"""
ADMET toxicity prediction using ADMET-AI (Chemprop-RDKit ensemble, 41 TDC datasets).
Falls back to RDKit heuristics if ADMET-AI is not installed.

Relaxed filter policy (v2):
  - hERG and Hepatotoxicity are the ONLY hard-fail properties.
    These are the two endpoints with the strongest clinical evidence for withdrawal/rejection.
  - Caco2 and Oral Bioavailability are INFORMATIONAL ONLY (hard_fail=False).
    Bioavailability_Ma is a binary classifier with ~30% FPR on diverse scaffolds —
    not reliable enough to discard otherwise good candidates.
  - BBB remains informational (unchanged).
  - "high" severity thresholds raised to match real-world approved-drug distributions:
    hERG >0.85 and hepatotox >0.85 are genuinely dangerous; 0.75 was over-aggressive.
"""

from typing import List, Tuple, Dict, Any
from backend.models.schemas import ADMETResult, ADMETFlagDetail


# ── ADMET-AI import with graceful fallback ────────────────────────────────────
try:
    from admet_ai import ADMETModel
    _ADMET_MODEL     = ADMETModel()
    _ADMET_AVAILABLE = True
except ImportError:
    _ADMET_AVAILABLE = False
    _ADMET_MODEL     = None


# ── Threshold config ──────────────────────────────────────────────────────────
#
# hard_fail=True  → compound is REMOVED from pipeline if severity=="high"
# hard_fail=False → flag is shown in report but compound continues downstream
#
# Severity escalation:
#   "moderate" = flag triggered but value is in a caution zone
#   "high"     = value is in the genuinely dangerous / very poor zone

THRESHOLDS: Dict[str, Any] = {
    "herg_inhibition": {
        "cutoff":         0.50,
        "direction":      "above",
        "threshold_str":  ">0.50",
        # Raised from 0.75 → 0.85.
        # Most approved oral drugs sit 0.50–0.75; only >0.85 signals true liability.
        "severity_high":  0.85,
        "implication":    "cardiac arrhythmia / QT prolongation risk",
        "recommendation": (
            "Redesign to remove hERG-liable substructure "
            "(basic amine + aromatic hydrophobic bulk). "
            "Consider bioisostere substitution or pKa reduction."
        ),
        "hard_fail":      True,   # ✅ keep — strongest clinical evidence
    },
    "hepatotoxicity": {
        "cutoff":         0.50,
        "direction":      "above",
        "threshold_str":  ">0.50",
        # Raised from 0.75 → 0.85 — same rationale as hERG.
        "severity_high":  0.85,
        "implication":    "liver toxicity risk — potential DILI (drug-induced liver injury)",
        "recommendation": (
            "Check for reactive metabolite formation. "
            "Consider metabolic soft-spot analysis; reduce lipophilicity (logP)."
        ),
        "hard_fail":      True,   # ✅ keep — #1 FDA rejection reason
    },
    "caco2_permeability": {
        "cutoff":         -5.15,
        "direction":      "below",
        "threshold_str":  "<-5.15 log cm/s",
        # Tightened: only truly impermeable compounds (< -6.5) get "high"
        "severity_high":  -6.5,
        "implication":    "poor intestinal absorption — low oral bioavailability predicted",
        "recommendation": (
            "Reduce polar surface area (PSA < 140 Å²) or HBD count. "
            "Consider prodrug strategy or permeation enhancer formulation."
        ),
        # ❌ Changed to False — Caco2_Wang has wide uncertainty intervals;
        # many drugs with poor Caco2 (e.g. bisphosphonates, aminoglycosides) are
        # clinically effective via non-Caco2 absorption routes.
        "hard_fail":      False,
    },
    "oral_bioavailability": {
        "cutoff":         0.30,
        "direction":      "below",
        "threshold_str":  "<0.30",
        # Only flag as "high" if essentially zero predicted bioavailability
        "severity_high":  0.10,
        "implication":    "poor systemic exposure after oral dosing",
        "recommendation": (
            "Evaluate first-pass metabolism and P-gp efflux. "
            "Consider IV or alternative delivery route, or structural optimisation."
        ),
        # ❌ Changed to False — Bioavailability_Ma (TDC) is a binary classifier
        # trained on a small dataset (~1k compounds); ~30% FPR on diverse scaffolds.
        # Useful as a flag, not reliable enough to hard-fail.
        "hard_fail":      False,
    },
    "bbb_penetration": {
        "cutoff":         0.30,
        "direction":      "below",
        "threshold_str":  "<0.30",
        "severity_high":  None,   # always "moderate"
        "implication":    "limited CNS penetration — may not reach brain target",
        "recommendation": (
            "Only relevant if CNS activity is required. "
            "Increase lipophilicity and reduce HBD count to improve CNS exposure."
        ),
        "hard_fail":      False,  # unchanged — informational only
    },
}


# ADMET-AI DataFrame column name → our schema field name
_COL_MAP = {
    "herg_inhibition":      "hERG",
    "caco2_permeability":   "Caco2_Wang",
    "bbb_penetration":      "BBB_Martins",
    "hepatotoxicity":       "DILI",
    "oral_bioavailability": "Bioavailability_Ma",
}

_PROP_LABELS = {
    "herg_inhibition":      "hERG Inhibition",
    "caco2_permeability":   "Caco-2 Permeability",
    "bbb_penetration":      "BBB Penetration",
    "hepatotoxicity":       "Hepatotoxicity",
    "oral_bioavailability": "Oral Bioavailability",
}

# Reverse map for _should_pass
_LABEL_TO_KEY = {v: k for k, v in _PROP_LABELS.items()}


# ── Public API ────────────────────────────────────────────────────────────────

def run_admet_batch(
    smiles_list: List[str],
) -> Tuple[List[str], List[ADMETResult]]:
    if not _ADMET_AVAILABLE:
        return _fallback_admet(smiles_list)
    try:
        df = _ADMET_MODEL.predict(smiles=smiles_list)
        return _parse_admet_dataframe(df, smiles_list)
    except Exception as e:
        print(f"[ADMET] ADMET-AI inference failed: {e}. Using fallback.")
        return _fallback_admet(smiles_list)


# ── Flag builder ──────────────────────────────────────────────────────────────

def _build_flags(
    herg:     float,
    caco2:    float,
    bbb:      float,
    hepato:   float,
    bioavail: float,
) -> Tuple[List[ADMETFlagDetail], List[str]]:
    values = {
        "herg_inhibition":      herg,
        "caco2_permeability":   caco2,
        "bbb_penetration":      bbb,
        "hepatotoxicity":       hepato,
        "oral_bioavailability": bioavail,
    }

    flags:   List[ADMETFlagDetail] = []
    summary: List[str]             = []

    for key, cfg in THRESHOLDS.items():
        val       = values[key]
        cutoff    = cfg["cutoff"]
        direction = cfg["direction"]

        triggered = (
            (direction == "above" and val > cutoff) or
            (direction == "below" and val < cutoff)
        )
        if not triggered:
            continue

        sev_thresh = cfg["severity_high"]
        if sev_thresh is None:
            severity = "moderate"
        elif direction == "above":
            severity = "high" if val > sev_thresh else "moderate"
        else:
            severity = "high" if val < sev_thresh else "moderate"

        # Non-hard-fail properties are capped at "moderate" regardless of value.
        # This prevents them from ever triggering _should_pass's hard-fail check.
        if not cfg["hard_fail"]:
            severity = "moderate"

        flag = ADMETFlagDetail(
            property_name  = _PROP_LABELS[key],
            value          = round(val, 3),
            threshold      = cfg["threshold_str"],
            direction      = direction,
            severity       = severity,
            implication    = cfg["implication"],
            recommendation = cfg["recommendation"],
        )
        flags.append(flag)
        summary.append(
            f"{_PROP_LABELS[key]}: {val:.2f} ({direction} {cfg['threshold_str']}) "
            f"[{severity}] — {cfg['implication']}"
        )

    return flags, summary


def _should_pass(flags: List[ADMETFlagDetail]) -> bool:
    """
    Fail only on HIGH-severity flags from hard_fail=True properties.
    Since non-hard-fail properties are capped at "moderate" in _build_flags,
    this check is effectively: hERG > 0.85 OR hepatotox > 0.85.
    """
    for f in flags:
        key = _LABEL_TO_KEY.get(f.property_name)
        if key and THRESHOLDS[key]["hard_fail"] and f.severity == "high":
            return False
    return True


# ── DataFrame parser ──────────────────────────────────────────────────────────

def _parse_admet_dataframe(
    df,
    smiles_list: List[str],
) -> Tuple[List[str], List[ADMETResult]]:
    passed_smiles: List[str]        = []
    all_results:   List[ADMETResult] = []

    for i, smi in enumerate(smiles_list):
        try:
            row = df.iloc[i]

            herg   = _safe_get(row, _COL_MAP["herg_inhibition"],      0.5)
            caco2  = _safe_get(row, _COL_MAP["caco2_permeability"],   -4.5)
            bbb    = _safe_get(row, _COL_MAP["bbb_penetration"],       0.5)
            hepato = _safe_get(row, _COL_MAP["hepatotoxicity"],        0.5)
            bioav  = _safe_get(row, _COL_MAP["oral_bioavailability"],  0.5)

            flags, summary = _build_flags(herg, caco2, bbb, hepato, bioav)
            passed         = _should_pass(flags)

            result = ADMETResult(
                passed               = passed,
                herg_inhibition      = round(herg,   4),
                caco2_permeability   = round(caco2,  4),
                bbb_penetration      = round(bbb,    4),
                hepatotoxicity       = round(hepato, 4),
                oral_bioavailability = round(bioav,  4),
                flags                = flags,
                flag_summary         = summary,
            )
            all_results.append(result)
            if passed:
                passed_smiles.append(smi)

        except Exception as e:
            print(f"[ADMET] Error parsing row {i} for SMILES '{smi}': {e}")
            all_results.append(_make_error_result(f"ADMET parse error: {e}"))

    return passed_smiles, all_results


# ── Helpers ───────────────────────────────────────────────────────────────────

def _safe_get(row, col_name: str, default: float) -> float:
    try:
        val = row.get(col_name, default)
        return float(val) if val is not None else default
    except (TypeError, ValueError):
        return default


def _make_error_result(reason: str) -> ADMETResult:
    err_flag = ADMETFlagDetail(
        property_name  = "Pipeline Error",
        value          = 0.0,
        threshold      = "N/A",
        direction      = "above",
        severity       = "high",
        implication    = reason,
        recommendation = "Re-submit the job or check the input SMILES.",
    )
    return ADMETResult(
        passed               = False,
        herg_inhibition      = 1.0,
        caco2_permeability   = -6.0,
        bbb_penetration      = 0.0,
        hepatotoxicity       = 1.0,
        oral_bioavailability = 0.0,
        flags                = [err_flag],
        flag_summary         = [f"Pipeline Error — {reason}"],
    )


# ── Fallback (no ADMET-AI installed) ─────────────────────────────────────────

def _fallback_admet(
    smiles_list: List[str],
) -> Tuple[List[str], List[ADMETResult]]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    results: List[ADMETResult] = []
    passed:  List[str]         = []

    warn_flag = ADMETFlagDetail(
        property_name  = "ADMET-AI Unavailable",
        value          = 0.0,
        threshold      = "N/A",
        direction      = "above",
        severity       = "moderate",
        implication    = "ADMET-AI not installed — values are rough RDKit heuristics only",
        recommendation = "Install admet-ai (pip install admet-ai) for accurate predictions.",
    )

    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                results.append(_make_error_result("Invalid SMILES in fallback"))
                continue

            logp       = Descriptors.MolLogP(mol)
            herg_rough = round(min(0.1 + max(0.0, (logp - 2.0) * 0.1), 1.0), 4)

            result = ADMETResult(
                passed               = True,
                herg_inhibition      = herg_rough,
                caco2_permeability   = -4.5,
                bbb_penetration      = 0.5,
                hepatotoxicity       = 0.3,
                oral_bioavailability = 0.6,
                flags                = [warn_flag],
                flag_summary         = [
                    "ADMET-AI Unavailable: 0.00 — rough heuristic estimates only [moderate]"
                ],
            )
            results.append(result)
            passed.append(smi)

        except Exception as e:
            results.append(_make_error_result(f"Fallback ADMET error: {e}"))

    return passed, results
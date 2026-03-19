from typing import List, Optional
from backend.models.schemas import CompoundResult


# ── Scoring weights (must sum to 100) ───────────────────────────────────────
# Each weight reflects how much that property matters for drug candidacy.
# Docking affinity is the primary driver; synthesis feasibility matters for
# real-world viability; ADMET safety is a strong secondary filter.

WEIGHTS = {
    "docking_affinity":    35,   # GNINA best affinity (lower kcal/mol = better)
    "cnn_score":           15,   # GNINA CNN pose quality (higher = better)
    "admet_safety":        20,   # Composite ADMET score (fewer flags = better)
    "solubility":          10,   # ESOL LogS solubility class
    "binding_prefilter":   10,   # DeepChem GNN predicted affinity
    "synthesis_ease":      10,   # Inverse of complexity_score (simpler = better)
}

assert sum(WEIGHTS.values()) == 100, "Scoring weights must sum to 100"

# Reference values for normalization
AFFINITY_BEST  = -12.0   # kcal/mol — excellent binder
AFFINITY_WORST =  -4.0   # kcal/mol — weak/no binding (fail threshold)

COMPLEXITY_BEST  =  2.0  # very simple synthesis
COMPLEXITY_WORST = 30.0  # extremely complex / infeasible

SOLUBILITY_SCORES = {
    "High":      1.0,
    "Medium":    0.7,
    "Low":       0.3,
    "Very Low":  0.0,
    "Invalid":   0.0,
}


def score_and_rank(compounds: List[CompoundResult]) -> List[CompoundResult]:
    """
    Computes a weighted composite score (0–100) for each compound
    and returns the list sorted by final_score descending (best first).

    Compounds missing a required score component (e.g. docking not run)
    receive 0 for that component — they sink to the bottom naturally.
    """
    for compound in compounds:
        compound.final_score = _compute_score(compound)

    ranked = sorted(compounds, key=lambda c: c.final_score or 0.0, reverse=True)
    for rank, compound in enumerate(ranked, start=1):
        compound.rank = rank

    return ranked


def _compute_score(c: CompoundResult) -> float:
    """
    Computes the composite score for a single CompoundResult.
    Each component is normalized to [0, 1] then multiplied by its weight.
    """
    total = 0.0

    # ── 1. Docking affinity (35 pts) ────────────────────────────────────────
    if c.docking and c.docking.best_affinity_kcal != 0.0:
        aff = c.docking.best_affinity_kcal
        norm = _normalize(aff, AFFINITY_WORST, AFFINITY_BEST, invert=False)
        total += norm * WEIGHTS["docking_affinity"]

    # ── 2. CNN score (15 pts) ────────────────────────────────────────────────
    if c.docking and c.docking.cnn_score > 0:
        # CNN score is already 0–1 range in GNINA
        norm = min(max(c.docking.cnn_score, 0.0), 1.0)
        total += norm * WEIGHTS["cnn_score"]

    # ── 3. ADMET safety (20 pts) ─────────────────────────────────────────────
    if c.admet:
        num_flags = len(c.admet.flags)
        # 0 flags = 1.0, 1 flag = 0.6, 2 flags = 0.3, 3+ flags = 0.0
        safety_scores = {0: 1.0, 1: 0.6, 2: 0.3}
        norm = safety_scores.get(num_flags, 0.0)
        # Also factor in hERG specifically (highest clinical risk)
        herg_penalty = max(0.0, c.admet.herg_inhibition - 0.3) * 2.0
        norm = max(0.0, norm - herg_penalty)
        total += norm * WEIGHTS["admet_safety"]

    # ── 4. Solubility (10 pts) ───────────────────────────────────────────────
    if c.lipinski:
        sol_norm = SOLUBILITY_SCORES.get(c.lipinski.solubility_class, 0.0)
        total += sol_norm * WEIGHTS["solubility"]

    # ── 5. Binding prefilter GNN score (10 pts) ──────────────────────────────
    if c.binding_prefilter:
        pre = c.binding_prefilter.predicted_affinity_kcal
        norm = _normalize(pre, AFFINITY_WORST, AFFINITY_BEST, invert=False)
        confidence_factor = c.binding_prefilter.confidence
        total += norm * confidence_factor * WEIGHTS["binding_prefilter"]

    # ── 6. Synthesis ease (10 pts) ───────────────────────────────────────────
    if c.retrosynthesis and c.retrosynthesis.feasible:
        comp = c.retrosynthesis.complexity_score
        norm = _normalize(comp, COMPLEXITY_WORST, COMPLEXITY_BEST, invert=True)
        total += norm * WEIGHTS["synthesis_ease"]
    elif c.retrosynthesis and not c.retrosynthesis.feasible:
        total += 0.0  # synthesis infeasible = 0 pts

    return round(min(max(total, 0.0), 100.0), 2)


def _normalize(
    value: float,
    worst: float,
    best: float,
    invert: bool = False,
) -> float:
    """
    Linearly normalizes a value to [0, 1] given worst and best reference points.
    invert=True means lower raw value = higher normalized score (e.g. complexity).
    """
    if best == worst:
        return 0.0
    norm = (value - worst) / (best - worst)
    if invert:
        norm = 1.0 - norm
    return min(max(norm, 0.0), 1.0)


def get_score_breakdown(compound: CompoundResult) -> dict:
    """
    Returns a human-readable breakdown of how each component contributed
    to the final score. Useful for the frontend expandable row.
    """
    c = compound
    breakdown = {}

    if c.docking and c.docking.best_affinity_kcal != 0.0:
        aff = c.docking.best_affinity_kcal
        norm = _normalize(aff, AFFINITY_WORST, AFFINITY_BEST)
        breakdown["docking_affinity"] = {
            "raw": f"{aff} kcal/mol",
            "contribution": round(norm * WEIGHTS["docking_affinity"], 2),
            "max_possible": WEIGHTS["docking_affinity"],
        }

    if c.docking and c.docking.cnn_score > 0:
        norm = min(max(c.docking.cnn_score, 0.0), 1.0)
        breakdown["cnn_score"] = {
            "raw": round(c.docking.cnn_score, 4),
            "contribution": round(norm * WEIGHTS["cnn_score"], 2),
            "max_possible": WEIGHTS["cnn_score"],
        }

    if c.admet:
        num_flags = len(c.admet.flags)
        safety_scores = {0: 1.0, 1: 0.6, 2: 0.3}
        norm = safety_scores.get(num_flags, 0.0)
        herg_penalty = max(0.0, c.admet.herg_inhibition - 0.3) * 2.0
        norm = max(0.0, norm - herg_penalty)
        breakdown["admet_safety"] = {
            "raw": f"{num_flags} flag(s): {', '.join(c.admet.flags) or 'None'}",
            "contribution": round(norm * WEIGHTS["admet_safety"], 2),
            "max_possible": WEIGHTS["admet_safety"],
        }

    if c.lipinski:
        sol_norm = SOLUBILITY_SCORES.get(c.lipinski.solubility_class, 0.0)
        breakdown["solubility"] = {
            "raw": f"{c.lipinski.solubility_class} (LogS={c.lipinski.logs})",
            "contribution": round(sol_norm * WEIGHTS["solubility"], 2),
            "max_possible": WEIGHTS["solubility"],
        }

    if c.retrosynthesis:
        if c.retrosynthesis.feasible:
            comp = c.retrosynthesis.complexity_score
            norm = _normalize(comp, COMPLEXITY_WORST, COMPLEXITY_BEST, invert=True)
            breakdown["synthesis_ease"] = {
                "raw": f"Complexity={comp}, Steps={c.retrosynthesis.num_steps}",
                "contribution": round(norm * WEIGHTS["synthesis_ease"], 2),
                "max_possible": WEIGHTS["synthesis_ease"],
            }
        else:
            breakdown["synthesis_ease"] = {
                "raw": "Synthesis infeasible",
                "contribution": 0.0,
                "max_possible": WEIGHTS["synthesis_ease"],
            }

    breakdown["final_score"] = c.final_score
    return breakdown

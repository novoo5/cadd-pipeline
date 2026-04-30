from typing import List, Optional
from backend.models.schemas import CompoundResult, PipelineSteps


# ── Scoring weights (must sum to 100) ─────────────────────────────────────────
WEIGHTS = {
    "docking_affinity":  50,
    "admet_safety":      20,
    "solubility":        10,
    "binding_prefilter": 10,
    "synthesis_ease":    10,
}

assert sum(WEIGHTS.values()) == 100, "Scoring weights must sum to 100"

AFFINITY_BEST  = -10.0
AFFINITY_WORST =  -3.0

COMPLEXITY_BEST  = 20.0
COMPLEXITY_WORST = 60.0

SOLUBILITY_SCORES = {
    "High":     1.0,
    "Medium":   0.7,
    "Low":      0.3,
    "Very Low": 0.0,
    "Invalid":  0.0,
}

SOLUBILITY_MULTIPLIERS = {
    "High":     1.00,
    "Medium":   1.00,
    "Low":      0.80,
    "Very Low": 0.55,
    "Invalid":  0.55,
}

MW_MIN = 200.0


def score_and_rank(
    compounds:      List[CompoundResult],
    pipeline_steps: Optional[PipelineSteps] = None,
) -> List[CompoundResult]:
    for compound in compounds:
        compound.final_score = _compute_score(compound, pipeline_steps)

    ranked = sorted(compounds, key=lambda c: c.final_score or 0.0, reverse=True)
    for rank, compound in enumerate(ranked, start=1):
        compound.rank = rank

    return ranked


def _docking_was_queued(c: CompoundResult) -> bool:
    """
    True if this compound made it through the binding prefilter but has no docking
    result — meaning docking was attempted but didn't complete (timeout/crash).
    Used only in get_score_breakdown() for human-readable explanation.
    """
    return (
        c.binding_prefilter is not None
        and c.binding_prefilter.passed
        and c.docking is None
    )


def _compute_score(
    c:              CompoundResult,
    pipeline_steps: Optional[PipelineSteps] = None,
) -> float:

    # Hard zero: docking was enabled in the pipeline but this compound has no result.
    # Covers timeouts, prefilter wipeouts, crashes, etc.
    # A compound that never docked cannot outrank compounds that did.
    if pipeline_steps is not None and pipeline_steps.docking and c.docking is None:
        return 0.0

    earned   = 0.0
    possible = 0.0

    # ── 1. Docking affinity (50 pts) ──────────────────────────────────────────
    if c.docking and c.docking.best_affinity_kcal != 0.0:
        aff      = c.docking.best_affinity_kcal
        earned  += _normalize(aff, AFFINITY_WORST, AFFINITY_BEST) * WEIGHTS["docking_affinity"]
        possible += WEIGHTS["docking_affinity"]
    elif _docking_was_queued(c):
        # Docking was attempted but didn't finish — award 0 pts but keep weight
        # in possible so the score isn't inflated by the missing component.
        possible += WEIGHTS["docking_affinity"]

    # ── 2. ADMET safety (20 pts) ──────────────────────────────────────────────
    if c.admet:
        num_flags    = len(c.admet.flags)
        safety_map   = {0: 1.0, 1: 0.6, 2: 0.3}
        norm         = safety_map.get(num_flags, 0.0)
        herg_penalty = max(0.0, c.admet.herg_inhibition - 0.3) * 2.0
        norm         = max(0.0, norm - herg_penalty)
        earned      += norm * WEIGHTS["admet_safety"]
        possible    += WEIGHTS["admet_safety"]

    # ── 3. Solubility (10 pts) ────────────────────────────────────────────────
    if c.lipinski:
        sol_norm  = SOLUBILITY_SCORES.get(c.lipinski.solubility_class, 0.0)
        earned   += sol_norm * WEIGHTS["solubility"]
        possible += WEIGHTS["solubility"]

    # ── 4. Binding prefilter GNN score (10 pts) ───────────────────────────────
    if c.binding_prefilter and c.binding_prefilter.confidence >= 0.4:
        pre      = c.binding_prefilter.predicted_affinity_kcal
        norm     = _normalize(pre, AFFINITY_WORST, AFFINITY_BEST)
        earned  += norm * WEIGHTS["binding_prefilter"]
        possible += WEIGHTS["binding_prefilter"]

    # ── 5. Synthesis ease (10 pts) ────────────────────────────────────────────
    if c.retrosynthesis and c.retrosynthesis.feasible:
        comp     = c.retrosynthesis.complexity_score
        norm     = _normalize(comp, COMPLEXITY_WORST, COMPLEXITY_BEST)
        earned  += norm * WEIGHTS["synthesis_ease"]
        possible += WEIGHTS["synthesis_ease"]
    elif c.retrosynthesis and not c.retrosynthesis.feasible:
        possible += WEIGHTS["synthesis_ease"]

    if possible == 0.0:
        return 0.0

    scaled = (earned / possible) * 100.0

    # ── MW fragment penalty ───────────────────────────────────────────────────
    if c.lipinski and c.lipinski.mw < MW_MIN:
        scaled *= 0.5

    # ── Solubility multiplier penalty ─────────────────────────────────────────
    if c.lipinski:
        sol_mult = SOLUBILITY_MULTIPLIERS.get(c.lipinski.solubility_class, 1.0)
        scaled  *= sol_mult

    return round(min(max(scaled, 0.0), 100.0), 2)


def _normalize(value: float, worst: float, best: float) -> float:
    if best == worst:
        return 0.0
    return min(max((value - worst) / (best - worst), 0.0), 1.0)


def get_score_breakdown(compound: CompoundResult) -> dict:
    c         = compound
    breakdown = {}

    # ── Docking ───────────────────────────────────────────────────────────────
    if c.docking and c.docking.best_affinity_kcal != 0.0:
        aff  = c.docking.best_affinity_kcal
        norm = _normalize(aff, AFFINITY_WORST, AFFINITY_BEST)
        breakdown["docking_affinity"] = {
            "raw":          f"{aff} kcal/mol",
            "contribution": round(norm * WEIGHTS["docking_affinity"], 2),
            "max_possible": WEIGHTS["docking_affinity"],
        }
    elif _docking_was_queued(c):
        breakdown["docking_affinity"] = {
            "raw":          "Docking did not complete (job timeout) — 0 pts, compound excluded from ranking",
            "contribution": 0.0,
            "max_possible": WEIGHTS["docking_affinity"],
        }
    elif c.docking is None:
        breakdown["docking_affinity"] = {
            "raw":          "No docking result — compound scored 0 and excluded from ranking",
            "contribution": 0.0,
            "max_possible": WEIGHTS["docking_affinity"],
        }

    # ── ADMET ─────────────────────────────────────────────────────────────────
    if c.admet:
        num_flags    = len(c.admet.flags)
        safety_map   = {0: 1.0, 1: 0.6, 2: 0.3}
        norm         = safety_map.get(num_flags, 0.0)
        herg_penalty = max(0.0, c.admet.herg_inhibition - 0.3) * 2.0
        norm         = max(0.0, norm - herg_penalty)
        flag_text    = ", ".join(c.admet.flag_summary) if c.admet.flag_summary else "None"
        breakdown["admet_safety"] = {
            "raw":          f"{num_flags} flag(s): {flag_text}",
            "contribution": round(norm * WEIGHTS["admet_safety"], 2),
            "max_possible": WEIGHTS["admet_safety"],
        }

    # ── Solubility ────────────────────────────────────────────────────────────
    if c.lipinski:
        sol_class = c.lipinski.solubility_class
        sol_norm  = SOLUBILITY_SCORES.get(sol_class, 0.0)
        sol_mult  = SOLUBILITY_MULTIPLIERS.get(sol_class, 1.0)
        mult_note = f" — {int((1 - sol_mult) * 100)}% final score penalty" if sol_mult < 1.0 else ""
        breakdown["solubility"] = {
            "raw":          f"{sol_class} (LogS={c.lipinski.logs:.2f}){mult_note}",
            "contribution": round(sol_norm * WEIGHTS["solubility"], 2),
            "max_possible": WEIGHTS["solubility"],
        }

    # ── Binding prefilter ─────────────────────────────────────────────────────
    if c.binding_prefilter and c.binding_prefilter.confidence >= 0.4:
        pre  = c.binding_prefilter.predicted_affinity_kcal
        norm = _normalize(pre, AFFINITY_WORST, AFFINITY_BEST)
        breakdown["binding_prefilter"] = {
            "raw":          f"{pre} kcal/mol (confidence={c.binding_prefilter.confidence:.2f})",
            "contribution": round(norm * WEIGHTS["binding_prefilter"], 2),
            "max_possible": WEIGHTS["binding_prefilter"],
        }

    # ── Retrosynthesis ────────────────────────────────────────────────────────
    if c.retrosynthesis:
        if c.retrosynthesis.feasible:
            comp = c.retrosynthesis.complexity_score
            norm = _normalize(comp, COMPLEXITY_WORST, COMPLEXITY_BEST)
            breakdown["synthesis_ease"] = {
                "raw":          f"Complexity={comp:.1f}, Steps={c.retrosynthesis.num_steps}",
                "contribution": round(norm * WEIGHTS["synthesis_ease"], 2),
                "max_possible": WEIGHTS["synthesis_ease"],
            }
        else:
            breakdown["synthesis_ease"] = {
                "raw":          "Synthesis infeasible",
                "contribution": 0.0,
                "max_possible": WEIGHTS["synthesis_ease"],
            }

    # ── Penalties ─────────────────────────────────────────────────────────────
    breakdown["mw_fragment_penalty"] = bool(c.lipinski and c.lipinski.mw < MW_MIN)
    breakdown["final_score"]         = c.final_score
    return breakdown
from typing import List, Optional
from backend.models.schemas import CompoundResult


# ── Scoring weights (must sum to 100) ───────────────────────────────────────
WEIGHTS = {
    "docking_affinity":    50,
    "admet_safety":        20,
    "solubility":          10,
    "binding_prefilter":   10,
    "synthesis_ease":      10,
}

assert sum(WEIGHTS.values()) == 100, "Scoring weights must sum to 100"

# ── Affinity range (realistic for Vina on drug-sized compounds) ──────────────
# FIX: was -12/-4 — too wide, making real scores compress into bottom 25%
AFFINITY_BEST  = -10.0   # excellent binding
AFFINITY_WORST =  -3.0   # poor binding (anything weaker = 0 pts)

# ── Complexity range (matches actual SA-derived values ~20–60) ───────────────
# FIX: was 2/30 — actual data shows 37–42, so everything was clamping to 10pts
COMPLEXITY_BEST  = 20.0   # simple to synthesise
COMPLEXITY_WORST = 60.0   # very hard to synthesise

# ── Solubility scores ────────────────────────────────────────────────────────
SOLUBILITY_SCORES = {
    "High":      1.0,
    "Medium":    0.7,
    "Low":       0.3,
    "Very Low":  0.0,
    "Invalid":   0.0,
}

# ── Minimum MW guard — penalise tiny fragments that pass Lipinski trivially ──
MW_MIN = 200.0   # anything below gets a 0.5× score multiplier


def score_and_rank(compounds: List[CompoundResult]) -> List[CompoundResult]:
    for compound in compounds:
        compound.final_score = _compute_score(compound)

    ranked = sorted(compounds, key=lambda c: c.final_score or 0.0, reverse=True)
    for rank, compound in enumerate(ranked, start=1):
        compound.rank = rank

    return ranked


def _compute_score(c: CompoundResult) -> float:
    total = 0.0

    # ── 1. Docking affinity (50 pts) ─────────────────────────────────────────
    if c.docking and c.docking.best_affinity_kcal != 0.0:
        aff = c.docking.best_affinity_kcal
        norm = _normalize(aff, AFFINITY_WORST, AFFINITY_BEST)
        total += norm * WEIGHTS["docking_affinity"]

    # ── 2. ADMET safety (20 pts) ──────────────────────────────────────────────
    if c.admet:
        num_flags = len(c.admet.flags)
        safety_scores = {0: 1.0, 1: 0.6, 2: 0.3}
        norm = safety_scores.get(num_flags, 0.0)
        # hERG penalty: only kicks in above 0.3 threshold
        herg_penalty = max(0.0, c.admet.herg_inhibition - 0.3) * 2.0
        norm = max(0.0, norm - herg_penalty)
        total += norm * WEIGHTS["admet_safety"]

    # ── 3. Solubility (10 pts) ────────────────────────────────────────────────
    if c.lipinski:
        sol_norm = SOLUBILITY_SCORES.get(c.lipinski.solubility_class, 0.0)
        total += sol_norm * WEIGHTS["solubility"]

    # ── 4. Binding prefilter GNN score (10 pts) ───────────────────────────────
    # FIX: don't multiply by confidence — it artificially caps this at ~5–8pts.
    # Instead, only use prefilter if confidence >= 0.4; treat it as 0 otherwise.
    if c.binding_prefilter and c.binding_prefilter.confidence >= 0.4:
        pre = c.binding_prefilter.predicted_affinity_kcal
        norm = _normalize(pre, AFFINITY_WORST, AFFINITY_BEST)
        total += norm * WEIGHTS["binding_prefilter"]

    # ── 5. Synthesis ease (10 pts) ────────────────────────────────────────────
    # FIX: removed invert=True — direction is handled by worst > best
    if c.retrosynthesis and c.retrosynthesis.feasible:
        comp = c.retrosynthesis.complexity_score
        norm = _normalize(comp, COMPLEXITY_WORST, COMPLEXITY_BEST)
        total += norm * WEIGHTS["synthesis_ease"]
    # infeasible or missing → 0 pts (no else needed)

    # ── MW fragment penalty ───────────────────────────────────────────────────
    # Tiny fragments (MW < 200) trivially pass Lipinski + ADMET — penalise them
    if c.lipinski and c.lipinski.molecular_weight < MW_MIN:
        total *= 0.5

    return round(min(max(total, 0.0), 100.0), 2)


def _normalize(value: float, worst: float, best: float) -> float:
    """Normalise value to [0, 1]. worst → 0.0, best → 1.0."""
    if best == worst:
        return 0.0
    norm = (value - worst) / (best - worst)
    return min(max(norm, 0.0), 1.0)


def get_score_breakdown(compound: CompoundResult) -> dict:
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

    if c.binding_prefilter and c.binding_prefilter.confidence >= 0.4:
        pre = c.binding_prefilter.predicted_affinity_kcal
        norm = _normalize(pre, AFFINITY_WORST, AFFINITY_BEST)
        breakdown["binding_prefilter"] = {
            "raw": f"{pre} kcal/mol (confidence={c.binding_prefilter.confidence:.2f})",
            "contribution": round(norm * WEIGHTS["binding_prefilter"], 2),
            "max_possible": WEIGHTS["binding_prefilter"],
        }

    if c.retrosynthesis:
        if c.retrosynthesis.feasible:
            comp = c.retrosynthesis.complexity_score
            norm = _normalize(comp, COMPLEXITY_WORST, COMPLEXITY_BEST)
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

    mw_penalty_applied = (
        c.lipinski and c.lipinski.molecular_weight < MW_MIN
    )
    breakdown["mw_fragment_penalty"] = mw_penalty_applied
    breakdown["final_score"] = c.final_score
    return breakdown

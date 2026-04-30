"""
Dynamic compound scoring and ranking.

v3 changes over v2:
  - Explicit weight redistribution via _resolve_active_weights():
    disabled pipeline steps lose their weight entirely; remaining
    components are scaled proportionally so effective weights always
    sum to 100. No denominator tricks — it's a first-class operation.
  - Kinder logS scoring: logS=-3 is genuinely soluble and should not
    score 3/10. New breakpoints: >=−4 = good, >=−5 = marginal, <−5 = bad.
  - _compute_score() uses resolved weights directly, not base weights.
  - Breakdown now shows both base_weight and effective_weight per component.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from backend.models.schemas import CompoundResult, PipelineSteps


# ── Scoring config ────────────────────────────────────────────────────────────


@dataclass
class ScoringConfig:
    """
    All scoring parameters in one place.

    Base weights define relative importance when ALL steps run.
    When a step is disabled or skipped, its weight is redistributed
    proportionally to the remaining active components via
    _resolve_active_weights(). The final scoring denominator is always 100.
    """

    # Base weights — must sum to 100
    weights: Dict[str, float] = field(
        default_factory=lambda: {
            "docking_affinity": 45,
            "admet_safety":     25,
            "drug_likeness":    10,
            "binding_prefilter": 12,
            "synthesis_ease":    8,
        }
    )

    # Maps component name → PipelineSteps attribute name
    # Used by _resolve_active_weights() to check which steps are enabled.
    component_to_step: Dict[str, str] = field(
        default_factory=lambda: {
            "docking_affinity":  "docking",
            "admet_safety":      "admet",
            "drug_likeness":     "drug_likeness",
            "binding_prefilter": "binding_prefilter",
            "synthesis_ease":    "retrosynthesis",
        }
    )

    # Docking / prefilter affinity normalization (kcal/mol)
    affinity_best:  float = -12.0
    affinity_worst: float = -3.0

    # Retrosynthesis complexity normalization
    complexity_best:  float = 15.0
    complexity_worst: float = 65.0

    # Minimum binding prefilter confidence to include in scoring at all
    prefilter_min_confidence: float = 0.25

    # MW floor below which final score is halved
    mw_fragment_floor: float = 200.0

    # ── ADMET ──
    admet_risk_scores: Dict[str, float] = field(
        default_factory=lambda: {
            "low":      1.00,
            "moderate": 0.55,
            "high":     0.00,
        }
    )
    admet_flag_penalty: Dict[str, float] = field(
        default_factory=lambda: {
            "high":     0.40,
            "moderate": 0.12,
        }
    )
    admet_hard_fail_multiplier: float = 2.0

    # ── Drug-likeness ──
    violation_scores: Dict[int, float] = field(
        default_factory=lambda: {
            0: 1.00,
            1: 0.80,
            2: 0.50,
            3: 0.20,
        }
    )

    # LogS breakpoints — (threshold, score): applies when logS >= threshold.
    # Sorted from least-negative to most-negative (highest to lowest logS).
    #
    # Scientific rationale:
    #   logS > -3  = soluble in any practical sense (~1 mg/mL+) → should score well.
    #   logS -3→-4 = still decent for oral delivery; many approved drugs sit here.
    #   logS -4→-5 = borderline; formulation may help but it's not a blocker.
    #   logS < -5  = genuinely poor solubility; starts mattering for bioavailability.
    #   logS < -6  = problematic; will likely need a formulation strategy.
    #
    # OLD scores gave logS=-3 → 0.30 (solubility class "Low" = 3/10).
    # That was way too harsh; the new curve gives logS=-3 → 0.90.
    logs_breakpoints: List[Tuple[float, float]] = field(
        default_factory=lambda: [
            (-1.0, 1.00),   # logS ≥ -1:  very soluble → perfect score
            (-3.0, 0.90),   # logS ≥ -3:  soluble → excellent (not 3/10 anymore)
            (-4.0, 0.75),   # logS ≥ -4:  decent → good ("below 4 is good")
            (-5.0, 0.45),   # logS ≥ -5:  borderline → mediocre ("5 starts being bad")
            (-6.0, 0.20),   # logS ≥ -6:  poor → low
        ]
    )
    logs_default_score: float = 0.05   # logS < -6: essentially insoluble

    # logP optimality windows
    logp_ideal_min:    float = 1.0
    logp_ideal_max:    float = 3.0
    logp_ok_min:       float = 0.0
    logp_ok_max:       float = 5.0
    logp_penalty_rate: float = 0.20   # score loss per unit above logp_ok_max

    # Drug-likeness sub-weights (must sum to 1.0)
    dl_violations_weight: float = 0.50
    dl_logs_weight:       float = 0.30
    dl_logp_weight:       float = 0.20

    def __post_init__(self) -> None:
        total = sum(self.weights.values())
        assert abs(total - 100) < 1e-6, \
            f"ScoringConfig.weights must sum to 100, got {total}"
        dl_sum = self.dl_violations_weight + self.dl_logs_weight + self.dl_logp_weight
        assert abs(dl_sum - 1.0) < 1e-6, \
            f"Drug-likeness sub-weights must sum to 1.0, got {dl_sum}"


DEFAULT_CONFIG = ScoringConfig()

_HARD_FAIL_ENDPOINTS = {"hERG Inhibition", "Hepatotoxicity"}


# ── Public API ────────────────────────────────────────────────────────────────


def score_and_rank(
    compounds:      List[CompoundResult],
    pipeline_steps: Optional[PipelineSteps] = None,
    config:         Optional[ScoringConfig] = None,
) -> List[CompoundResult]:
    cfg = config or DEFAULT_CONFIG
    effective_weights = _resolve_active_weights(pipeline_steps, cfg)

    for compound in compounds:
        compound.final_score = _compute_score(compound, pipeline_steps, cfg, effective_weights)

    ranked = sorted(compounds, key=lambda c: c.final_score or 0.0, reverse=True)
    for rank, compound in enumerate(ranked, start=1):
        compound.rank = rank

    return ranked


def get_score_breakdown(
    compound:       CompoundResult,
    pipeline_steps: Optional[PipelineSteps] = None,
    config:         Optional[ScoringConfig] = None,
) -> dict:
    cfg = config or DEFAULT_CONFIG
    effective_weights = _resolve_active_weights(pipeline_steps, cfg)
    return _build_breakdown(compound, cfg, effective_weights)


def get_effective_weights(
    pipeline_steps: Optional[PipelineSteps] = None,
    config:         Optional[ScoringConfig] = None,
) -> Dict[str, float]:
    """
    Returns the redistributed weights for the given pipeline configuration.
    Call this from the frontend to show users which components are active
    and what weight each carries after redistribution.
    """
    cfg = config or DEFAULT_CONFIG
    return _resolve_active_weights(pipeline_steps, cfg)


# ── Weight redistribution ─────────────────────────────────────────────────────


def _resolve_active_weights(
    pipeline_steps: Optional[PipelineSteps],
    cfg: ScoringConfig,
) -> Dict[str, float]:
    """
    Compute per-component effective weights after redistributing the weight
    of any disabled pipeline step to the remaining active components.

    Example: if docking (45 pts) is disabled:
      - active base total = 55 pts
      - each active weight is scaled by 100/55
      - docking_affinity effective weight = 0.0

    The returned dict always sums to 100.0 across active components.
    """
    if pipeline_steps is None:
        # No pipeline context → use base weights as-is
        return {k: float(w) for k, w in cfg.weights.items()}

    active: Dict[str, float] = {}
    disabled: List[str] = []

    for component, step_attr in cfg.component_to_step.items():
        is_enabled = getattr(pipeline_steps, step_attr, True)
        if is_enabled:
            active[component] = cfg.weights[component]
        else:
            disabled.append(component)

    if not active:
        return {k: 0.0 for k in cfg.weights}

    active_total = sum(active.values())
    scale = 100.0 / active_total

    result: Dict[str, float] = {}
    for component in cfg.weights:
        if component in active:
            result[component] = round(active[component] * scale, 4)
        else:
            result[component] = 0.0

    return result


# ── Core scoring ──────────────────────────────────────────────────────────────


def _compute_score(
    c:                CompoundResult,
    pipeline_steps:   Optional[PipelineSteps],
    cfg:              ScoringConfig,
    effective_weights: Dict[str, float],
) -> float:
    """
    Hard zero: docking was enabled but no result arrived.
    A compound that never completed docking cannot outrank one that did.
    """
    if pipeline_steps is not None and pipeline_steps.docking and c.docking is None:
        return 0.0

    earned = 0.0
    denominator = 0.0

    # Each scorer returns (score_0_to_1) or None if data is absent.
    # We only include a component in the denominator when:
    #   (a) its step is enabled (effective_weight > 0), AND
    #   (b) data for it actually exists on this compound.
    # If a step is enabled but data is absent (e.g. prefilter low confidence),
    # the weight is still included in the denominator to prevent inflation.

    def _add(key: str, score: Optional[float], *, data_present: bool) -> None:
        nonlocal earned, denominator
        w = effective_weights.get(key, 0.0)
        if w == 0.0:
            return  # step disabled; skip entirely
        denominator += w
        if score is not None and data_present:
            earned += score * w

    _add("docking_affinity",  _score_docking(c, cfg),            data_present=c.docking is not None)
    _add("admet_safety",      _score_admet(c, cfg),              data_present=c.admet is not None)
    _add("drug_likeness",     _score_drug_likeness(c, cfg),      data_present=c.lipinski is not None)
    _add("binding_prefilter", _score_binding_prefilter(c, cfg),  data_present=c.binding_prefilter is not None)
    _add("synthesis_ease",    _score_synthesis(c, cfg),          data_present=c.retrosynthesis is not None)

    if denominator == 0.0:
        return 0.0

    scaled = (earned / denominator) * 100.0

    # MW fragment penalty
    if c.lipinski and c.lipinski.mw < cfg.mw_fragment_floor:
        scaled *= 0.5

    return round(min(max(scaled, 0.0), 100.0), 2)


# ── Component scorers ─────────────────────────────────────────────────────────


def _score_docking(c: CompoundResult, cfg: ScoringConfig) -> Optional[float]:
    if c.docking is None or c.docking.best_affinity_kcal == 0.0:
        return None
    return _normalize(c.docking.best_affinity_kcal, cfg.affinity_worst, cfg.affinity_best)


def _score_admet(c: CompoundResult, cfg: ScoringConfig) -> Optional[float]:
    if c.admet is None:
        return None

    risk = getattr(c.admet, "overall_risk", "moderate")
    base = cfg.admet_risk_scores.get(risk, 0.5)

    if base == 0.0:
        return 0.0

    cumulative_penalty = 0.0
    for flag in c.admet.flags:
        pen = cfg.admet_flag_penalty.get(flag.severity, 0.0)
        if flag.property_name in _HARD_FAIL_ENDPOINTS:
            pen *= cfg.admet_hard_fail_multiplier
        cumulative_penalty += pen

    return max(0.0, base - min(cumulative_penalty, base))


def _score_drug_likeness(c: CompoundResult, cfg: ScoringConfig) -> Optional[float]:
    if c.lipinski is None:
        return None

    # 1. Lipinski violations sub-score
    violations  = _count_lipinski_violations(c.lipinski)
    viol_score  = cfg.violation_scores.get(violations, 0.0)

    # 2. logS continuous sub-score
    logs       = c.lipinski.logs
    logs_score = cfg.logs_default_score
    for threshold, score in sorted(cfg.logs_breakpoints, reverse=True):
        if logs >= threshold:
            logs_score = score
            break

    # 3. logP optimality sub-score
    logp = c.lipinski.logp
    if cfg.logp_ideal_min <= logp <= cfg.logp_ideal_max:
        logp_score = 1.00
    elif cfg.logp_ok_min <= logp <= cfg.logp_ok_max:
        logp_score = 0.75
    elif logp > cfg.logp_ok_max:
        logp_score = max(0.0, 1.0 - (logp - cfg.logp_ok_max) * cfg.logp_penalty_rate)
    else:
        logp_score = max(0.50, 1.0 + logp * 0.10)

    combined = (
        viol_score * cfg.dl_violations_weight
        + logs_score * cfg.dl_logs_weight
        + logp_score * cfg.dl_logp_weight
    )
    return min(max(combined, 0.0), 1.0)


def _score_binding_prefilter(c: CompoundResult, cfg: ScoringConfig) -> Optional[float]:
    if c.binding_prefilter is None:
        return None
    if c.binding_prefilter.confidence < cfg.prefilter_min_confidence:
        return None

    base = _normalize(
        c.binding_prefilter.predicted_affinity_kcal, cfg.affinity_worst, cfg.affinity_best
    )

    conf     = c.binding_prefilter.confidence
    min_conf = cfg.prefilter_min_confidence
    if min_conf >= 1.0:
        conf_weight = 1.0
    else:
        conf_weight = 0.40 + 0.60 * (conf - min_conf) / (1.0 - min_conf)

    return min(1.0, base * conf_weight)


def _score_synthesis(c: CompoundResult, cfg: ScoringConfig) -> Optional[float]:
    if c.retrosynthesis is None:
        return None
    if not c.retrosynthesis.feasible:
        return 0.0

    complexity_score = _normalize(
        c.retrosynthesis.complexity_score, cfg.complexity_worst, cfg.complexity_best
    )

    steps = c.retrosynthesis.num_steps
    if steps <= 2:
        step_score = 1.00
    elif steps <= 3:
        step_score = 0.90
    elif steps <= 5:
        step_score = 0.70
    elif steps <= 7:
        step_score = 0.45
    else:
        step_score = 0.20

    return complexity_score * 0.60 + step_score * 0.40


# ── Breakdown ─────────────────────────────────────────────────────────────────


def _build_breakdown(
    c:                 CompoundResult,
    cfg:               ScoringConfig,
    effective_weights: Dict[str, float],
) -> dict:
    bd: dict = {}

    def _entry(
        key:          str,
        score:        Optional[float],
        raw:          str,
        explanation:  str,
        data_present: bool = True,
    ) -> dict:
        ew = effective_weights.get(key, 0.0)
        bw = cfg.weights.get(key, 0.0)
        s  = score if (score is not None and data_present) else 0.0
        return {
            "score_0_to_1":    round(s, 3),
            "raw":             raw,
            "contribution":    round(s * ew, 2),
            "base_weight":     bw,
            "effective_weight": round(ew, 2),
            "explanation":     explanation,
            "step_active":     ew > 0.0,
        }

    # ── Docking ───────────────────────────────────────────────────────────────
    dock_score = _score_docking(c, cfg)
    if c.docking is not None:
        if dock_score is not None:
            aff = c.docking.best_affinity_kcal
            bd["docking_affinity"] = _entry(
                "docking_affinity",
                dock_score,
                f"{aff:.2f} kcal/mol",
                (
                    f"Affinity {aff:.2f} normalised between "
                    f"worst={cfg.affinity_worst} and best={cfg.affinity_best} kcal/mol."
                ),
            )
        else:
            bd["docking_affinity"] = _entry(
                "docking_affinity",
                0.0,
                "Affinity = 0.0 — docking error.",
                "Zero affinity value excluded from scoring.",
            )
    elif _docking_was_attempted(c):
        bd["docking_affinity"] = _entry(
            "docking_affinity",
            0.0,
            "Passed prefilter — docking did not complete (timeout/crash).",
            "Compound scored 0 and excluded from ranking.",
        )

    # ── ADMET ─────────────────────────────────────────────────────────────────
    admet_score = _score_admet(c, cfg)
    if c.admet is not None:
        risk       = getattr(c.admet, "overall_risk", "moderate")
        flags_text = "; ".join(c.admet.flag_summary[:3]) if c.admet.flag_summary else "None"
        bd["admet_safety"] = _entry(
            "admet_safety",
            admet_score,
            f"Risk: {risk} | {len(c.admet.flags)} flag(s) | {flags_text}",
            (
                f"Base={cfg.admet_risk_scores.get(risk, 0):.2f} (risk={risk}) minus "
                f"per-flag penalties (moderate={cfg.admet_flag_penalty['moderate']}, "
                f"high={cfg.admet_flag_penalty['high']}; "
                f"×{cfg.admet_hard_fail_multiplier} for hERG/hepato)."
            ),
        )

    # ── Drug-likeness ─────────────────────────────────────────────────────────
    dl_score = _score_drug_likeness(c, cfg)
    if c.lipinski is not None:
        violations = _count_lipinski_violations(c.lipinski)
        bd["drug_likeness"] = _entry(
            "drug_likeness",
            dl_score,
            (
                f"Violations={violations} | "
                f"LogS={c.lipinski.logs:.2f} ({c.lipinski.solubility_class}) | "
                f"LogP={c.lipinski.logp:.2f} | MW={c.lipinski.mw:.1f} Da"
            ),
            (
                f"Lipinski violations ({cfg.dl_violations_weight*100:.0f}%) + "
                f"logS ({cfg.dl_logs_weight*100:.0f}%) + "
                f"logP optimality ({cfg.dl_logp_weight*100:.0f}%). "
                "LogS ≥ -4 = good; ≥ -5 = marginal; < -5 = poor."
            ),
        )

    # ── Binding prefilter ─────────────────────────────────────────────────────
    pre_score = _score_binding_prefilter(c, cfg)
    if c.binding_prefilter is not None:
        conf    = c.binding_prefilter.confidence
        pre_aff = c.binding_prefilter.predicted_affinity_kcal
        bd["binding_prefilter"] = _entry(
            "binding_prefilter",
            pre_score,
            f"{pre_aff:.2f} kcal/mol (confidence={conf:.2f})",
            (
                "Normalised affinity scaled by confidence weight. "
                f"Min confidence to include: {cfg.prefilter_min_confidence}. "
                "Low confidence → partial contribution rather than hard exclusion."
            ),
            data_present=(pre_score is not None),
        )

    # ── Synthesis ease ────────────────────────────────────────────────────────
    synth_score = _score_synthesis(c, cfg)
    if c.retrosynthesis is not None:
        bd["synthesis_ease"] = _entry(
            "synthesis_ease",
            synth_score,
            (
                f"Steps={c.retrosynthesis.num_steps} | "
                f"Complexity={c.retrosynthesis.complexity_score:.1f} | "
                f"Feasible={c.retrosynthesis.feasible}"
            ),
            "60% normalised complexity + 40% step-count penalty.",
        )

    # ── Weight redistribution info ────────────────────────────────────────────
    bd["weight_redistribution"] = {
        "base_weights": {k: cfg.weights[k] for k in cfg.weights},
        "effective_weights": {k: round(v, 2) for k, v in effective_weights.items()},
        "disabled_components": [
            k for k, v in effective_weights.items() if v == 0.0
        ],
        "note": (
            "Effective weights always sum to 100. "
            "Disabled steps have their base weight redistributed "
            "proportionally to active components."
        ),
    }

    # ── Penalties ─────────────────────────────────────────────────────────────
    mw = c.lipinski.mw if c.lipinski else None
    bd["penalties"] = {
        "mw_fragment_penalty_applied": bool(mw is not None and mw < cfg.mw_fragment_floor),
        "mw_value":    mw,
        "mw_floor":    cfg.mw_fragment_floor,
        "effect":      "50% final score halved" if (mw is not None and mw < cfg.mw_fragment_floor) else "none",
    }

    bd["final_score"] = c.final_score
    return bd


# ── Helpers ───────────────────────────────────────────────────────────────────


def _normalize(value: float, worst: float, best: float) -> float:
    if best == worst:
        return 0.0
    return min(max((value - worst) / (best - worst), 0.0), 1.0)


def _count_lipinski_violations(lipinski) -> int:
    violations = 0
    if lipinski.mw   > 500:  violations += 1
    if lipinski.logp > 5:    violations += 1
    if lipinski.hbd  > 5:    violations += 1
    if lipinski.hba  > 10:   violations += 1
    return violations


def _docking_was_attempted(c: CompoundResult) -> bool:
    return (
        c.binding_prefilter is not None
        and c.binding_prefilter.passed
        and c.docking is None
    )
"""
ADMET toxicity prediction using ADMET-AI (Chemprop-RDKit ensemble, 41 TDC datasets).
Falls back to RDKit heuristics if ADMET-AI is not installed.

Updated version:
 - Keeps current 5 decision-driving endpoints.
 - Adds configurable threshold support via JobRequest.admet_config.
 - Adds ADMET narrative / write-up blocks for frontend display.
 - Exposes a generic endpoint table + extra ADMET-AI properties for richer UI.
 - Remains backward-compatible with old callers: run_admet_batch(smiles_list) still works.
"""

from typing import List, Tuple, Dict, Any, Optional

from backend.models.schemas import (
    ADMETResult,
    ADMETFlagDetail,
    ADMETEndpointValue,
    ADMETNarrativeBlock,
    ADMETTuningConfig,
    ADMETThresholdConfig,
)

# ── ADMET-AI import with graceful fallback ────────────────────────────────────
try:
    from admet_ai import ADMETModel
    _ADMET_MODEL = ADMETModel()
    _ADMET_AVAILABLE = True
except ImportError:
    _ADMET_AVAILABLE = False
    _ADMET_MODEL = None


# ── Core endpoint metadata ────────────────────────────────────────────────────

_ENDPOINT_META: Dict[str, Dict[str, Any]] = {
    "herg_inhibition": {
        "column": "hERG",
        "label": "hERG Inhibition",
        "category": "toxicity",
        "default": 0.5,
        "unit": "probability",
    },
    "caco2_permeability": {
        "column": "Caco2_Wang",
        "label": "Caco-2 Permeability",
        "category": "absorption",
        "default": -4.5,
        "unit": "log cm/s",
    },
    "bbb_penetration": {
        "column": "BBB_Martins",
        "label": "BBB Penetration",
        "category": "distribution",
        "default": 0.5,
        "unit": "probability",
    },
    "hepatotoxicity": {
        "column": "DILI",
        "label": "Hepatotoxicity",
        "category": "toxicity",
        "default": 0.5,
        "unit": "probability",
    },
    "oral_bioavailability": {
        "column": "Bioavailability_Ma",
        "label": "Oral Bioavailability",
        "category": "absorption",
        "default": 0.5,
        "unit": "probability",
    },
}

_PROP_LABELS = {k: v["label"] for k, v in _ENDPOINT_META.items()}
_LABEL_TO_KEY = {v["label"]: k for k, v in _ENDPOINT_META.items()}
_CORE_COLUMNS = {v["column"] for v in _ENDPOINT_META.values()}


# ── Default thresholds ────────────────────────────────────────────────────────

DEFAULT_THRESHOLDS: Dict[str, ADMETThresholdConfig] = {
    "herg_inhibition": ADMETThresholdConfig(
        cutoff=0.50,
        direction="above",
        severity_high=0.85,
        hard_fail=True,
        threshold_str=">0.50",
        implication="cardiac arrhythmia / QT prolongation risk",
        recommendation=(
            "Redesign to remove hERG-liable substructure "
            "(basic amine + aromatic hydrophobic bulk). "
            "Consider bioisostere substitution or pKa reduction."
        ),
    ),
    "hepatotoxicity": ADMETThresholdConfig(
        cutoff=0.50,
        direction="above",
        severity_high=0.85,
        hard_fail=True,
        threshold_str=">0.50",
        implication="liver toxicity risk — potential DILI (drug-induced liver injury)",
        recommendation=(
            "Check for reactive metabolite formation. "
            "Consider metabolic soft-spot analysis; reduce lipophilicity (logP)."
        ),
    ),
    "caco2_permeability": ADMETThresholdConfig(
        cutoff=-5.15,
        direction="below",
        severity_high=-6.5,
        hard_fail=False,
        threshold_str="<-5.15 log cm/s",
        implication="poor intestinal absorption — low oral bioavailability predicted",
        recommendation=(
            "Reduce polar surface area (PSA < 140 Å²) or HBD count. "
            "Consider prodrug strategy or permeation enhancer formulation."
        ),
    ),
    "oral_bioavailability": ADMETThresholdConfig(
        cutoff=0.30,
        direction="below",
        severity_high=0.10,
        hard_fail=False,
        threshold_str="<0.30",
        implication="poor systemic exposure after oral dosing",
        recommendation=(
            "Evaluate first-pass metabolism and P-gp efflux. "
            "Consider IV or alternative delivery route, or structural optimisation."
        ),
    ),
    "bbb_penetration": ADMETThresholdConfig(
        cutoff=0.30,
        direction="below",
        severity_high=None,
        hard_fail=False,
        threshold_str="<0.30",
        implication="limited CNS penetration — may not reach brain target",
        recommendation=(
            "Only relevant if CNS activity is required. "
            "Increase lipophilicity and reduce HBD count to improve CNS exposure."
        ),
    ),
}

# Optional preset patches. Manual UI edits should set preset="custom".
PRESET_PATCHES: Dict[str, Dict[str, Dict[str, Any]]] = {
    "balanced": {},
    "custom": {},
    "cns": {
        "bbb_penetration": {
            "cutoff": 0.60,
            "severity_high": 0.30,
            "threshold_str": "<0.60",
        }
    },
    "oral": {
        "caco2_permeability": {
            "cutoff": -5.00,
            "severity_high": -6.00,
            "threshold_str": "<-5.00 log cm/s",
        },
        "oral_bioavailability": {
            "cutoff": 0.40,
            "severity_high": 0.20,
            "threshold_str": "<0.40",
        },
    },
}


# ── Public API ────────────────────────────────────────────────────────────────

def run_admet_batch(
    smiles_list: List[str],
    admet_config: Optional[ADMETTuningConfig] = None,
) -> Tuple[List[str], List[ADMETResult]]:
    thresholds = _merge_thresholds(admet_config)

    if not _ADMET_AVAILABLE:
        return _fallback_admet(smiles_list, thresholds)

    try:
        df = _ADMET_MODEL.predict(smiles=smiles_list)
        return _parse_admet_dataframe(df, smiles_list, thresholds)
    except Exception as e:
        print(f"[ADMET] ADMET-AI inference failed: {e}. Using fallback.")
        return _fallback_admet(smiles_list, thresholds)


# ── Threshold resolution ──────────────────────────────────────────────────────

def _get_default_thresholds() -> Dict[str, ADMETThresholdConfig]:
    return {k: v.model_copy(deep=True) for k, v in DEFAULT_THRESHOLDS.items()}


def _models_differ(a: ADMETThresholdConfig, b: ADMETThresholdConfig) -> bool:
    return a.model_dump() != b.model_dump()


def _apply_preset(
    preset: str,
    thresholds: Dict[str, ADMETThresholdConfig],
) -> Dict[str, ADMETThresholdConfig]:
    patch = PRESET_PATCHES.get(preset, {})
    updated = {k: v.model_copy(deep=True) for k, v in thresholds.items()}

    for key, updates in patch.items():
        if key not in updated:
            continue
        merged = updated[key].model_dump()
        merged.update(updates)
        updated[key] = ADMETThresholdConfig(**merged)

    return updated


def _merge_thresholds(
    admet_config: Optional[ADMETTuningConfig],
) -> Dict[str, ADMETThresholdConfig]:
    thresholds = _get_default_thresholds()

    if admet_config is None:
        return thresholds

    thresholds = _apply_preset(admet_config.preset, thresholds)

    config_map = {
        "herg_inhibition": admet_config.herg_inhibition,
        "hepatotoxicity": admet_config.hepatotoxicity,
        "caco2_permeability": admet_config.caco2_permeability,
        "oral_bioavailability": admet_config.oral_bioavailability,
        "bbb_penetration": admet_config.bbb_penetration,
    }

    # If preset == balanced/custom, use supplied config directly.
    # If preset != balanced/custom, only overwrite preset values when config differs from default.
    for key, override in config_map.items():
        if admet_config.preset in {"balanced", "custom"}:
            thresholds[key] = override.model_copy(deep=True)
        else:
            if _models_differ(override, DEFAULT_THRESHOLDS[key]):
                thresholds[key] = override.model_copy(deep=True)

    return thresholds


# ── Flag builder ──────────────────────────────────────────────────────────────

def _format_threshold(cfg: ADMETThresholdConfig) -> str:
    if cfg.threshold_str:
        return cfg.threshold_str

    sign = ">" if cfg.direction == "above" else "<"
    return f"{sign}{cfg.cutoff:.2f}"


def _is_triggered(value: float, cfg: ADMETThresholdConfig) -> bool:
    if not cfg.enabled:
        return False
    if cfg.direction == "above":
        return value > cfg.cutoff
    return value < cfg.cutoff


def _compute_severity(value: float, cfg: ADMETThresholdConfig) -> str:
    if cfg.severity_high is None:
        return "moderate"

    if cfg.direction == "above":
        severity = "high" if value > cfg.severity_high else "moderate"
    else:
        severity = "high" if value < cfg.severity_high else "moderate"

    if not cfg.hard_fail:
        return "moderate"

    return severity


def _build_flags(
    values: Dict[str, float],
    thresholds: Dict[str, ADMETThresholdConfig],
) -> Tuple[List[ADMETFlagDetail], List[str]]:
    flags: List[ADMETFlagDetail] = []
    summary: List[str] = []

    for key, cfg in thresholds.items():
        val = values.get(key)
        if val is None or not cfg.enabled:
            continue

        triggered = _is_triggered(val, cfg)
        if not triggered:
            continue

        severity = _compute_severity(val, cfg)
        threshold_str = _format_threshold(cfg)

        implication = cfg.implication or "Potential ADMET liability."
        recommendation = cfg.recommendation or "Review this endpoint in context."

        flag = ADMETFlagDetail(
            property_name=_PROP_LABELS[key],
            value=round(val, 3),
            threshold=threshold_str,
            direction=cfg.direction,
            severity=severity,
            implication=implication,
            recommendation=recommendation,
        )
        flags.append(flag)
        summary.append(
            f"{_PROP_LABELS[key]}: {val:.2f} ({cfg.direction} {threshold_str}) "
            f"[{severity}] — {implication}"
        )

    return flags, summary


def _should_pass(
    flags: List[ADMETFlagDetail],
    thresholds: Dict[str, ADMETThresholdConfig],
) -> Tuple[bool, List[str]]:
    hard_fail_reasons: List[str] = []

    for f in flags:
        key = _LABEL_TO_KEY.get(f.property_name)
        if key and thresholds[key].hard_fail and f.severity == "high":
            hard_fail_reasons.append(
                f"{f.property_name} triggered hard fail at {f.value:.3f} "
                f"against threshold {f.threshold}"
            )

    if hard_fail_reasons:
        return False, hard_fail_reasons

    return True, ["No high-severity hard-fail ADMET liability detected."]


def _overall_risk_from_flags(flags: List[ADMETFlagDetail]) -> str:
    if any(f.severity == "high" for f in flags):
        return "high"
    if any(f.severity == "moderate" for f in flags):
        return "moderate"
    return "low"


# ── DataFrame parser ──────────────────────────────────────────────────────────

def _parse_admet_dataframe(
    df,
    smiles_list: List[str],
    thresholds: Dict[str, ADMETThresholdConfig],
) -> Tuple[List[str], List[ADMETResult]]:
    passed_smiles: List[str] = []
    all_results: List[ADMETResult] = []

    for i, smi in enumerate(smiles_list):
        try:
            row = df.iloc[i]

            core_values = {
                key: _safe_get(row, meta["column"], meta["default"])
                for key, meta in _ENDPOINT_META.items()
            }

            flags, summary = _build_flags(core_values, thresholds)
            passed, decision_basis = _should_pass(flags, thresholds)
            overall_risk = _overall_risk_from_flags(flags)

            extra_properties = _extract_extra_properties(row)
            endpoint_table = _build_endpoint_table(core_values, flags, thresholds, extra_properties)
            narrative = _build_narrative(core_values, flags, passed, overall_risk, extra_properties)

            result = ADMETResult(
                passed=passed,
                overall_risk=overall_risk,
                herg_inhibition=round(core_values["herg_inhibition"], 4),
                caco2_permeability=round(core_values["caco2_permeability"], 4),
                bbb_penetration=round(core_values["bbb_penetration"], 4),
                hepatotoxicity=round(core_values["hepatotoxicity"], 4),
                oral_bioavailability=round(core_values["oral_bioavailability"], 4),
                flags=flags,
                flag_summary=summary,
                decision_basis=decision_basis,
                narrative=narrative,
                endpoint_table=endpoint_table,
                extra_properties=extra_properties,
                thresholds_used={k: v.model_copy(deep=True) for k, v in thresholds.items()},
            )
            all_results.append(result)

            if passed:
                passed_smiles.append(smi)

        except Exception as e:
            print(f"[ADMET] Error parsing row {i} for SMILES '{smi}': {e}")
            all_results.append(_make_error_result(f"ADMET parse error: {e}", thresholds))

    return passed_smiles, all_results


# ── Endpoint table + narrative ────────────────────────────────────────────────

def _build_endpoint_table(
    core_values: Dict[str, float],
    flags: List[ADMETFlagDetail],
    thresholds: Dict[str, ADMETThresholdConfig],
    extra_properties: Dict[str, float],
) -> List[ADMETEndpointValue]:
    table: List[ADMETEndpointValue] = []
    flag_map = {f.property_name: f for f in flags}

    for key, value in core_values.items():
        meta = _ENDPOINT_META[key]
        label = meta["label"]
        flag = flag_map.get(label)
        cfg = thresholds[key]

        table.append(
            ADMETEndpointValue(
                key=key,
                label=label,
                category=meta["category"],
                value=round(value, 4),
                display_value=_format_value(value),
                unit=meta["unit"],
                interpretation=_interpret_core_endpoint(key, value),
                threshold_applied=_format_threshold(cfg),
                severity=flag.severity if flag else None,
                hard_fail=cfg.hard_fail,
                triggered=flag is not None,
            )
        )

    for col_name, value in sorted(extra_properties.items()):
        table.append(
            ADMETEndpointValue(
                key=col_name,
                label=col_name.replace("_", " "),
                category=_categorize_extra_column(col_name),
                value=round(value, 4),
                display_value=_format_value(value),
                unit=None,
                interpretation="Additional ADMET-AI endpoint returned by the model.",
                threshold_applied=None,
                severity=None,
                hard_fail=False,
                triggered=False,
            )
        )

    return table


def _build_narrative(
    core_values: Dict[str, float],
    flags: List[ADMETFlagDetail],
    passed: bool,
    overall_risk: str,
    extra_properties: Dict[str, float],
) -> List[ADMETNarrativeBlock]:
    blocks: List[ADMETNarrativeBlock] = []

    verdict = "passed" if passed else "failed"
    blocks.append(
        ADMETNarrativeBlock(
            title="Overall ADMET call",
            level="good" if passed else "danger",
            body=(
                f"This compound {verdict} the current ADMET decision rules with an overall "
                f"risk level of {overall_risk}. "
                f"hERG={core_values['herg_inhibition']:.3f}, "
                f"DILI={core_values['hepatotoxicity']:.3f}, "
                f"Caco-2={core_values['caco2_permeability']:.3f}, "
                f"Oral BA={core_values['oral_bioavailability']:.3f}, "
                f"BBB={core_values['bbb_penetration']:.3f}."
            ),
        )
    )

    if flags:
        top_lines = "; ".join(
            f"{f.property_name} [{f.severity}] — {f.implication}"
            for f in flags[:4]
        )
        blocks.append(
            ADMETNarrativeBlock(
                title="Key liabilities",
                level="warning" if passed else "danger",
                body=top_lines,
            )
        )
    else:
        blocks.append(
            ADMETNarrativeBlock(
                title="Key liabilities",
                level="good",
                body="No configured ADMET thresholds were triggered for this compound.",
            )
        )

    exposure_note = (
        f"Caco-2 permeability is {core_values['caco2_permeability']:.3f}, "
        f"predicted oral bioavailability is {core_values['oral_bioavailability']:.3f}, "
        f"and BBB penetration is {core_values['bbb_penetration']:.3f}. "
        f"These should be interpreted relative to the intended route and target tissue."
    )
    blocks.append(
        ADMETNarrativeBlock(
            title="Exposure and distribution",
            level="info",
            body=exposure_note,
        )
    )

    if extra_properties:
        blocks.append(
            ADMETNarrativeBlock(
                title="Extended ADMET profile",
                level="info",
                body=(
                    f"The model returned {len(extra_properties)} additional ADMET/PK properties "
                    f"beyond the 5 core decision endpoints. These can be shown in an expanded "
                    f"profile table in the frontend."
                ),
            )
        )

    return blocks


def _interpret_core_endpoint(key: str, value: float) -> str:
    if key == "herg_inhibition":
        if value > 0.85:
            return "High predicted cardiac liability."
        if value > 0.50:
            return "Moderate hERG concern; monitor carefully."
        return "No major hERG signal under current rule."

    if key == "hepatotoxicity":
        if value > 0.85:
            return "High predicted DILI/hepatotoxicity liability."
        if value > 0.50:
            return "Moderate liver-risk signal."
        return "No major hepatotoxicity signal under current rule."

    if key == "caco2_permeability":
        if value < -6.5:
            return "Very poor intestinal permeability."
        if value < -5.15:
            return "Suboptimal permeability for oral delivery."
        return "Permeability looks acceptable under current rule."

    if key == "oral_bioavailability":
        if value < 0.10:
            return "Very low predicted oral exposure."
        if value < 0.30:
            return "Suboptimal oral bioavailability."
        return "Oral bioavailability looks acceptable under current rule."

    if key == "bbb_penetration":
        if value < 0.30:
            return "Limited CNS penetration predicted."
        return "BBB penetration may be acceptable if CNS exposure is needed."

    return "No interpretation available."


def _categorize_extra_column(col_name: str) -> str:
    upper = col_name.upper()

    if any(x in upper for x in ["CACO2", "BIOAVAIL", "HIA", "PAMPA", "PGLYCOPROTEIN", "PGP"]):
        return "absorption"
    if any(x in upper for x in ["BBB", "PPBR", "VDSS", "VD", "FU"]):
        return "distribution"
    if any(x in upper for x in ["CYP", "MICROSOME", "UGT"]):
        return "metabolism"
    if any(x in upper for x in ["CLEARANCE", "HALF_LIFE", "T12", "T1_2"]):
        return "excretion"
    if any(x in upper for x in ["HERG", "DILI", "AMES", "TOX", "CARCIN", "SKIN", "RESPIR"]):
        return "toxicity"
    if any(x in upper for x in ["LOGP", "LOGS", "TPSA", "MW", "HBA", "HBD"]):
        return "physchem"

    return "other"


def _format_value(value: float) -> str:
    return f"{value:.3f}"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _safe_get(row, col_name: str, default: float) -> float:
    try:
        val = row.get(col_name, default)
        if val is None:
            return default
        out = float(val)
        if out != out:  # NaN check
            return default
        return out
    except (TypeError, ValueError):
        return default


def _extract_extra_properties(row) -> Dict[str, float]:
    props: Dict[str, float] = {}

    for col_name, val in row.items():
        if col_name in _CORE_COLUMNS:
            continue

        try:
            fval = float(val)
            if fval != fval:  # NaN
                continue
            props[str(col_name)] = round(fval, 4)
        except (TypeError, ValueError):
            continue

    return props


def _make_error_result(
    reason: str,
    thresholds: Optional[Dict[str, ADMETThresholdConfig]] = None,
) -> ADMETResult:
    err_flag = ADMETFlagDetail(
        property_name="Pipeline Error",
        value=0.0,
        threshold="N/A",
        direction="above",
        severity="high",
        implication=reason,
        recommendation="Re-submit the job or check the input SMILES.",
    )

    return ADMETResult(
        passed=False,
        overall_risk="high",
        herg_inhibition=1.0,
        caco2_permeability=-6.0,
        bbb_penetration=0.0,
        hepatotoxicity=1.0,
        oral_bioavailability=0.0,
        flags=[err_flag],
        flag_summary=[f"Pipeline Error — {reason}"],
        decision_basis=[reason],
        narrative=[
            ADMETNarrativeBlock(
                title="ADMET pipeline error",
                level="danger",
                body=reason,
            )
        ],
        endpoint_table=[],
        extra_properties={},
        thresholds_used=thresholds or {},
    )


# ── Fallback (no ADMET-AI installed) ─────────────────────────────────────────

def _fallback_admet(
    smiles_list: List[str],
    thresholds: Dict[str, ADMETThresholdConfig],
) -> Tuple[List[str], List[ADMETResult]]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    results: List[ADMETResult] = []
    passed: List[str] = []

    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                results.append(_make_error_result("Invalid SMILES in fallback", thresholds))
                continue

            logp = Descriptors.MolLogP(mol)
            herg_rough = round(min(0.1 + max(0.0, (logp - 2.0) * 0.1), 1.0), 4)

            core_values = {
                "herg_inhibition": herg_rough,
                "caco2_permeability": -4.5,
                "bbb_penetration": 0.5,
                "hepatotoxicity": 0.3,
                "oral_bioavailability": 0.6,
            }

            warn_flag = ADMETFlagDetail(
                property_name="ADMET-AI Unavailable",
                value=0.0,
                threshold="N/A",
                direction="above",
                severity="moderate",
                implication="ADMET-AI not installed — values are rough RDKit heuristics only",
                recommendation="Install admet-ai (pip install admet-ai) for accurate predictions.",
            )

            flags, summary = _build_flags(core_values, thresholds)
            flags = [warn_flag] + flags
            summary = ["ADMET-AI unavailable — rough heuristic estimates only."] + summary

            did_pass, decision_basis = _should_pass(flags, thresholds)
            overall_risk = _overall_risk_from_flags(flags)

            result = ADMETResult(
                passed=did_pass,
                overall_risk=overall_risk,
                herg_inhibition=core_values["herg_inhibition"],
                caco2_permeability=core_values["caco2_permeability"],
                bbb_penetration=core_values["bbb_penetration"],
                hepatotoxicity=core_values["hepatotoxicity"],
                oral_bioavailability=core_values["oral_bioavailability"],
                flags=flags,
                flag_summary=summary,
                decision_basis=decision_basis,
                narrative=[
                    ADMETNarrativeBlock(
                        title="Fallback ADMET mode",
                        level="warning",
                        body=(
                            "ADMET-AI was unavailable, so this profile uses rough RDKit-based "
                            "heuristics and default placeholder values for non-hERG endpoints."
                        ),
                    )
                ],
                endpoint_table=_build_endpoint_table(core_values, [], thresholds, {}),
                extra_properties={},
                thresholds_used={k: v.model_copy(deep=True) for k, v in thresholds.items()},
            )

            results.append(result)
            if did_pass:
                passed.append(smi)

        except Exception as e:
            results.append(_make_error_result(f"Fallback ADMET error: {e}", thresholds))

    return passed, results
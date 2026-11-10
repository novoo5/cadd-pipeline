"use client";

import { useState } from "react";
import {
    ChevronDown, ChevronUp, Info, Zap, FlaskConical,
    Lock, X, Clock, ShieldAlert, Settings2, Dna, Pill, Brain,
    Activity, Layers, Target, AlertTriangle, MapPin,
} from "lucide-react";
import type {
    PipelineSteps,
    DockingSpeed,
    BindingSiteMode,
    BindingSiteCoords,
    BindingSiteResidues,
    SolubilityFilterMode,
    ADMETTuningConfig,
    ADMETPreset,
} from "@/lib/api";
import { ADMET_PRESET_DEFAULTS } from "@/lib/api";

// ── Props ─────────────────────────────────────────────────────────────────────

interface AdvancedSettingsProps {
    numAnalogues: number;
    onNumAnaloguesChange: (v: number) => void;
    pipelineSteps: PipelineSteps;
    onPipelineStepsChange: (v: PipelineSteps) => void;
    dockingSpeed: DockingSpeed;
    onDockingSpeedChange: (v: DockingSpeed) => void;
    maxDockingCompounds: number;
    onMaxDockingCompoundsChange: (v: number) => void;
    bindingSiteMode: BindingSiteMode;
    onBindingSiteModeChange: (v: BindingSiteMode) => void;
    bindingSiteCoords: BindingSiteCoords;
    onBindingSiteCoordsChange: (v: BindingSiteCoords) => void;
    bindingSiteResidues: BindingSiteResidues;
    onBindingSiteResiduesChange: (v: BindingSiteResidues) => void;
    directScoreOnly: boolean;
    onDirectScoreOnlyChange: (v: boolean) => void;
    mwMin: number;
    mwMax: number;
    onMwMinChange: (v: number) => void;
    onMwMaxChange: (v: number) => void;
    maxLipinskiViolations: number | null;
    onMaxLipinskiViolationsChange: (v: number | null) => void;
    solubilityFilter: SolubilityFilterMode;
    onSolubilityFilterChange: (v: SolubilityFilterMode) => void;
    toxicityReportOnly: boolean;
    onToxicityReportOnlyChange: (v: boolean) => void;
    lockedScaffoldSmarts: string;
    onLockedScaffoldSmartsChange: (v: string) => void;
    admetConfig: ADMETTuningConfig;
    onAdmetConfigChange: (v: ADMETTuningConfig) => void;
}

// ── Static config ─────────────────────────────────────────────────────────────

const STEPS: { key: keyof PipelineSteps; label: string; desc: string; locked?: boolean; icon: React.ReactNode }[] = [
    { key: "drug_likeness", label: "Drug-likeness Filter", desc: "Lipinski RO5 + ESOL solubility (RDKit)", locked: true, icon: <Activity className="w-3.5 h-3.5" /> },
    { key: "admet", label: "ADMET Toxicity Filter", desc: "hERG, hepatotox, Caco-2, bioavailability (ADMET-AI)", icon: <ShieldAlert className="w-3.5 h-3.5" /> },
    { key: "binding_prefilter", label: "ML Binding Pre-filter", desc: "GNN affinity ranking — sends top N for docking (DeepChem AttentiveFP)", icon: <Layers className="w-3.5 h-3.5" /> },
    { key: "docking", label: "Molecular Docking", desc: "AutoDock Vina with CNN rescoring", icon: <Target className="w-3.5 h-3.5" /> },
    { key: "retrosynthesis", label: "Retrosynthesis", desc: "SA Score + feasibility (AiZynthFinder)", icon: <Dna className="w-3.5 h-3.5" /> },
];

const DOCKING_SPEEDS: { value: DockingSpeed; label: string; tag: string; detail: string; minPerCompound: number }[] = [
    { value: "fast", label: "Fast", tag: "=8", detail: "~5 min/cpd", minPerCompound: 5 },
    { value: "balanced", label: "Balanced", tag: "=16", detail: "~10 min/cpd", minPerCompound: 10 },
    { value: "thorough", label: "Thorough", tag: "=32", detail: "~20 min/cpd", minPerCompound: 20 },
];

const DOCKING_COUNT_PRESETS = [10, 25, 50] as const;

const VIOLATION_OPTIONS: { value: string; label: string; sub: string }[] = [
    { value: "1", label: "Strict", sub: "≤1 violation — classic RO5" },
    { value: "2", label: "Relaxed", sub: "≤2 violations" },
    { value: "3", label: "Lenient", sub: "≤3 violations" },
    { value: "null", label: "Disabled", sub: "No Lipinski filter" },
];

const SOLUBILITY_OPTIONS: { value: SolubilityFilterMode; label: string; sub: string }[] = [
    { value: "soluble_only", label: "Soluble only", sub: "logS > −3 · aqueous assay-ready" },
    { value: "allow_slightly", label: "Slightly soluble", sub: "logS > −5 · may need formulation" },
    { value: "all", label: "No filter", sub: "No solubility restriction" },
];

const SCAFFOLD_PRESETS = [
    { label: "Quinolone", smarts: "c1ccc2cc1C(=O)c1ccccc1N2", hint: "Fluoroquinolone bicyclic ring" },
    { label: "Benzimidazole", smarts: "c1ccc2nHcnc2c1", hint: "Fused benzimidazole core" },
    { label: "Flavone", smarts: "O=c1cc(-c2ccccc2)oc2ccccc12", hint: "2-phenyl-4H-chromen-4-one" },
    { label: "Purine", smarts: "c1nc2nHcnc2n1", hint: "Adenine-like, kinase hinge" },
    { label: "Indole", smarts: "c1ccc2nHccc2c1", hint: "Bicyclic indole ring system" },
];

const ANALOGUE_PRESETS = [25, 50, 100, 500, 1000] as const;

const ADMET_PRESETS: { value: ADMETPreset; label: string; desc: string; icon: React.ReactNode; accent: string }[] = [
    { value: "balanced", label: "Balanced", desc: "Standard early discovery thresholds", icon: <Settings2 className="w-4 h-4" />, accent: "teal" },
    { value: "oral", label: "Oral-focused", desc: "Stricter Caco-2 & bioavailability", icon: <Pill className="w-4 h-4" />, accent: "blue" },
    { value: "cns", label: "CNS-focused", desc: "Higher BBB penetration required", icon: <Brain className="w-4 h-4" />, accent: "violet" },
    { value: "custom", label: "Custom", desc: "Manual threshold adjustments", icon: <Dna className="w-4 h-4" />, accent: "amber" },
];

const PRESET_ACCENT: Record<string, string> = {
    teal: "border-teal-600/70   bg-teal-950/40   text-teal-300   [&_svg]:text-teal-400",
    blue: "border-blue-600/70   bg-blue-950/40   text-blue-300   [&_svg]:text-blue-400",
    violet: "border-violet-600/70 bg-violet-950/40 text-violet-300 [&_svg]:text-violet-400",
    amber: "border-amber-600/70  bg-amber-950/40  text-amber-300  [&_svg]:text-amber-400",
};

const BINDING_MODES: { value: BindingSiteMode; label: string; hint: string }[] = [
    { value: "auto", label: "Auto-detect", hint: "Detects site from co-crystallised ligand in the PDB file" },
    { value: "coordinates", label: "Manual Coords", hint: "Provide Cartesian centre + box size directly" },
    { value: "residues", label: "Residue Range", hint: "Define pocket by chain ID and residue number range" },
];

// ── ADMET endpoint config ─────────────────────────────────────────────────────

// Fix for TypeScript: CSS custom properties need this union type.
// Without it, `"--thumb-border": color` causes "jsx element can't have
// multiple attributes with the same name" in strict mode.
type CSSWithVars = React.CSSProperties & {
    [K in `--${string}`]?: string | number;
};

type ADMETEndpointKey =
    | "herg_inhibition"
    | "hepatotoxicity"
    | "caco2_permeability"
    | "oral_bioavailability"
    | "bbb_penetration";

interface ADMETEndpointCfg {
    key: ADMETEndpointKey;
    label: string;
    subtitle: string;
    safeNote: string;
    sliderField: "cutoff" | "severity_high";
    min: number;
    max: number;
    step: number;
    color: string;
    thumbColor: string;
    hardFail: boolean;
    failLabel: string;
    defaultGetter: (cfg: ADMETTuningConfig) => number;
    formatVal: (v: number) => string;
}

const ADMET_ENDPOINT_CONFIGS: ADMETEndpointCfg[] = [
    {
        key: "herg_inhibition",
        label: "hERG Inhibition",
        subtitle: "Cardiac arrhythmia · QT prolongation",
        safeNote: "Lower is safer",
        sliderField: "severity_high",
        min: 0.50, max: 0.95, step: 0.05,
        color: "rgb(248 113 113 / 0.7)",
        thumbColor: "#f87171",
        hardFail: true,
        failLabel: "Hard-fail if above",
        defaultGetter: (cfg) => cfg.herg_inhibition?.severity_high ?? 0.85,
        formatVal: (v) => v.toFixed(2),
    },
    {
        key: "hepatotoxicity",
        label: "Hepatotoxicity (DILI)",
        subtitle: "Drug-induced liver injury probability",
        safeNote: "Lower is safer",
        sliderField: "severity_high",
        min: 0.50, max: 0.95, step: 0.05,
        color: "rgb(248 113 113 / 0.7)",
        thumbColor: "#f87171",
        hardFail: true,
        failLabel: "Hard-fail if above",
        defaultGetter: (cfg) => cfg.hepatotoxicity?.severity_high ?? 0.85,
        formatVal: (v) => v.toFixed(2),
    },
    {
        key: "caco2_permeability",
        label: "Caco-2 Permeability",
        subtitle: "Intestinal absorption — log cm/s",
        safeNote: "Higher (less negative) is better",
        sliderField: "cutoff",
        min: -7.0, max: -4.0, step: 0.25,
        color: "rgb(96 165 250 / 0.7)",
        thumbColor: "#60a5fa",
        hardFail: false,
        failLabel: "Flag if below",
        defaultGetter: (cfg) => cfg.caco2_permeability?.cutoff ?? -5.15,
        formatVal: (v) => v.toFixed(2),
    },
    {
        key: "oral_bioavailability",
        label: "Oral Bioavailability",
        subtitle: "Fraction reaching systemic circulation",
        safeNote: "Higher is better",
        sliderField: "cutoff",
        min: 0.05, max: 0.60, step: 0.05,
        color: "rgb(96 165 250 / 0.7)",
        thumbColor: "#60a5fa",
        hardFail: false,
        failLabel: "Flag if below",
        defaultGetter: (cfg) => cfg.oral_bioavailability?.cutoff ?? 0.30,
        formatVal: (v) => v.toFixed(2),
    },
    {
        key: "bbb_penetration",
        label: "BBB Penetration",
        subtitle: "Blood-brain barrier probability",
        safeNote: "CNS drugs only — informational otherwise",
        sliderField: "cutoff",
        min: 0.05, max: 0.70, step: 0.05,
        color: "rgb(167 139 250 / 0.7)",
        thumbColor: "#a78bfa",
        hardFail: false,
        failLabel: "Flag if below",
        defaultGetter: (cfg) => cfg.bbb_penetration?.cutoff ?? 0.30,
        formatVal: (v) => v.toFixed(2),
    },
];

// ── Helpers ───────────────────────────────────────────────────────────────────

function formatDockingTime(compounds: number, minPerCompound: number): string {
    const total = compounds * minPerCompound;
    if (total < 60) return `~${total} min`;
    const h = Math.floor(total / 60);
    const m = total % 60;
    return m === 0 ? `~${h}h` : `~${h}h ${m}m`;
}

function dockingTimeSeverity(compounds: number, minPerCompound: number): "ok" | "warn" | "heavy" {
    const total = compounds * minPerCompound;
    if (total <= 60) return "ok";
    if (total <= 300) return "warn";
    return "heavy";
}

// Returns CSSWithVars so --thumb-border is a valid key — no "as" casting needed at call site
function sliderFillStyle(value: number, min: number, max: number, color: string, thumbColor: string): CSSWithVars {
    const pct = ((value - min) / (max - min)) * 100;
    return {
        background: `linear-gradient(to right, ${color} 0%, ${color} ${pct}%, rgb(55 65 81) ${pct}%, rgb(55 65 81) 100%)`,
        "--thumb-border": thumbColor,
    };
}

// ── Sub-components ────────────────────────────────────────────────────────────

function SectionHeader({
    icon,
    label,
    aside,
}: {
    icon: React.ReactNode;
    label: string;
    aside?: React.ReactNode;
}) {
    return (
        <div className="flex items-center justify-between mb-5">
            <div className="flex items-center gap-2.5">
                <span className="text-gray-500">{icon}</span>
                <span className="text-[11px] font-bold uppercase tracking-[0.12em] text-gray-400">
                    {label}
                </span>
            </div>
            {aside && <div className="flex items-center">{aside}</div>}
        </div>
    );
}

function SectionDivider() {
    return <div className="h-px bg-gradient-to-r from-transparent via-gray-700/60 to-transparent my-1" />;
}

function FieldLabel({ children, hint }: { children: React.ReactNode; hint?: string }) {
    return (
        <div className="flex items-center gap-1.5 mb-2.5">
            <span className="text-xs font-semibold text-gray-300">{children}</span>
            {hint && (
                <span className="group relative cursor-help">
                    <Info className="w-3 h-3 text-gray-600 group-hover:text-gray-400 transition-colors" />
                    <span className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 w-52 bg-gray-900 border border-gray-700 rounded-lg px-3 py-2 text-xs text-gray-300 invisible group-hover:visible z-50 pointer-events-none leading-relaxed shadow-xl">
                        {hint}
                    </span>
                </span>
            )}
        </div>
    );
}

function SegmentedPill<T extends string | number>({
    options,
    value,
    onChange,
    renderLabel,
}: {
    options: T[];
    value: T;
    onChange: (v: T) => void;
    renderLabel?: (v: T) => React.ReactNode;
}) {
    return (
        <div className="inline-flex items-center gap-1 bg-gray-900/60 border border-gray-700/60 rounded-xl p-1">
            {options.map((opt) => (
                <button
                    key={String(opt)}
                    type="button"
                    onClick={() => onChange(opt)}
                    className={`px-3.5 py-1.5 rounded-lg text-xs font-semibold transition-all ${value === opt
                        ? "bg-teal-600 text-white shadow-sm shadow-teal-900/50"
                        : "text-gray-400 hover:text-gray-200 hover:bg-gray-700/50"
                        }`}
                >
                    {renderLabel ? renderLabel(opt) : String(opt)}
                </button>
            ))}
        </div>
    );
}

function RunModeCard({
    checked,
    onChange,
    icon,
    title,
    desc,
    activeColor,
}: {
    checked: boolean;
    onChange: (v: boolean) => void;
    icon: React.ReactNode;
    title: string;
    desc: string;
    activeColor: string;
}) {
    return (
        <label className={`flex items-start gap-4 p-4 rounded-2xl border cursor-pointer transition-all select-none ${checked ? activeColor : "bg-gray-900/20 border-gray-800 hover:border-gray-600 hover:bg-gray-800/30"
            }`}>
            <input
                type="checkbox"
                checked={checked}
                onChange={(e) => onChange(e.target.checked)}
                className="sr-only"
            />
            <div className={`mt-0.5 w-5 h-5 rounded-md border-2 flex items-center justify-center shrink-0 transition-all ${checked ? "bg-current border-current" : "border-gray-600 bg-gray-800"
                }`}>
                {checked && (
                    <svg className="w-3 h-3 text-white" viewBox="0 0 12 12" fill="none">
                        <path d="M2 6l3 3 5-5" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
                    </svg>
                )}
            </div>
            <div className="min-w-0">
                <p className="text-sm font-semibold text-gray-100 flex items-center gap-2">
                    {icon} {title}
                </p>
                <p className="text-xs text-gray-500 mt-1.5 leading-relaxed">{desc}</p>
            </div>
        </label>
    );
}

// ── Component ─────────────────────────────────────────────────────────────────

export default function AdvancedSettings({
    numAnalogues, onNumAnaloguesChange,
    pipelineSteps, onPipelineStepsChange,
    dockingSpeed, onDockingSpeedChange,
    maxDockingCompounds, onMaxDockingCompoundsChange,
    bindingSiteMode, onBindingSiteModeChange,
    bindingSiteCoords, onBindingSiteCoordsChange,
    bindingSiteResidues, onBindingSiteResiduesChange,
    directScoreOnly, onDirectScoreOnlyChange,
    mwMin, mwMax, onMwMinChange, onMwMaxChange,
    maxLipinskiViolations, onMaxLipinskiViolationsChange,
    solubilityFilter, onSolubilityFilterChange,
    toxicityReportOnly, onToxicityReportOnlyChange,
    lockedScaffoldSmarts, onLockedScaffoldSmartsChange,
    admetConfig, onAdmetConfigChange,
}: AdvancedSettingsProps) {
    const [open, setOpen] = useState(false);
    const [scaffoldPresetOpen, setScaffoldPresetOpen] = useState(false);
    const [dockingCountCustom, setDockingCountCustom] = useState(false);

    const toxOnly = toxicityReportOnly && !directScoreOnly;
    const bothModes = toxicityReportOnly && directScoreOnly;
    const analoguesOff = directScoreOnly;
    const dockingDisabled = toxOnly || !pipelineSteps.docking;

    const toggleStep = (key: keyof PipelineSteps) => {
        if (key === "drug_likeness") return;
        onPipelineStepsChange({ ...pipelineSteps, [key]: !pipelineSteps[key] });
    };

    const violationVal = maxLipinskiViolations === null ? "null" : String(maxLipinskiViolations);
    const handleViolation = (raw: string) =>
        onMaxLipinskiViolationsChange(raw === "null" ? null : parseInt(raw, 10));

    const handleAdmetPreset = (preset: ADMETPreset) => {
        const defaults = ADMET_PRESET_DEFAULTS[preset];
        onAdmetConfigChange({ ...admetConfig, ...defaults, preset });
    };

    // Generic handler for all 5 ADMET endpoints
    const handleAdmetSlider = (
        endpoint: ADMETEndpointKey,
        field: "cutoff" | "severity_high",
        value: number,
    ) => {
        const current = (admetConfig[endpoint as keyof ADMETTuningConfig] as Record<string, number> | undefined)
            ?? (ADMET_PRESET_DEFAULTS["balanced"] as Record<string, Record<string, number>>)[endpoint]
            ?? {};
        onAdmetConfigChange({
            ...admetConfig,
            preset: "custom",
            [endpoint]: { ...current, [field]: value },
        });
    };

    // Toggle endpoint enabled/disabled
    const handleAdmetToggle = (endpoint: ADMETEndpointKey) => {
        const current = (admetConfig[endpoint as keyof ADMETTuningConfig] as Record<string, unknown> | undefined)
            ?? (ADMET_PRESET_DEFAULTS["balanced"] as Record<string, Record<string, unknown>>)[endpoint]
            ?? {};
        onAdmetConfigChange({
            ...admetConfig,
            preset: "custom",
            [endpoint]: { ...current, enabled: !(current["enabled"] ?? true) },
        });
    };

    const activeSpeed = DOCKING_SPEEDS.find(s => s.value === dockingSpeed)!;
    const dockTimeStr = formatDockingTime(maxDockingCompounds, activeSpeed.minPerCompound);
    const dockTimeSev = dockingTimeSeverity(maxDockingCompounds, activeSpeed.minPerCompound);
    const isCustomCount = !DOCKING_COUNT_PRESETS.includes(maxDockingCompounds as typeof DOCKING_COUNT_PRESETS[number]);

    const activeBadges = [
        directScoreOnly && { label: "Direct Score", color: "bg-violet-900/50 border-violet-700/50 text-violet-300" },
        toxicityReportOnly && { label: bothModes ? "Tox + Direct" : "Tox Only", color: "bg-red-900/50 border-red-700/50 text-red-300" },
        lockedScaffoldSmarts.trim() && { label: "Scaffold Lock", color: "bg-amber-900/50 border-amber-700/50 text-amber-300" },
        admetConfig.preset !== "balanced" && { label: `${admetConfig.preset} ADMET`, color: "bg-teal-900/50 border-teal-700/50 text-teal-300" },
        bindingSiteMode !== "auto" && { label: bindingSiteMode === "coordinates" ? "Manual Coords" : "Residue Range", color: "bg-sky-900/50 border-sky-700/50 text-sky-300" },
    ].filter(Boolean) as { label: string; color: string }[];

    // Hard-fail and advisory endpoint groups
    const hardFailEndpoints = ADMET_ENDPOINT_CONFIGS.filter(c => c.hardFail);
    const advisoryEndpoints = ADMET_ENDPOINT_CONFIGS.filter(c => !c.hardFail);

    return (
        <div className="card mt-4">
            {/* Header toggle */}
            <button
                type="button"
                onClick={() => setOpen(!open)}
                className="w-full flex items-center justify-between text-sm font-semibold text-gray-300 hover:text-white transition-colors group"
            >
                <span className="flex items-center gap-3 flex-wrap">
                    <Settings2 className="w-4 h-4 text-gray-500 group-hover:text-gray-400 transition-colors" />
                    <span>Advanced Settings</span>
                    {activeBadges.map(b => (
                        <span key={b.label} className={`inline-flex items-center px-2 py-0.5 rounded-full text-[11px] font-semibold border ${b.color}`}>
                            {b.label}
                        </span>
                    ))}
                </span>
                <span className="flex items-center gap-1.5 text-xs text-gray-600 group-hover:text-gray-400 transition-colors shrink-0">
                    {open ? <ChevronUp className="w-4 h-4" /> : <ChevronDown className="w-4 h-4" />}
                </span>
            </button>

            {open && (
                <div className="mt-8 space-y-10 animate-in fade-in slide-in-from-top-2 duration-200">

                    {/* SECTION 1 — Run Mode */}
                    <section>
                        <SectionHeader icon={<Zap className="w-4 h-4" />} label="Run Mode" />
                        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                            <RunModeCard
                                checked={directScoreOnly}
                                onChange={onDirectScoreOnlyChange}
                                icon={<Zap className="w-3.5 h-3.5 text-violet-400" />}
                                title="Direct Score Only"
                                desc="Skip analogue generation. Runs the input compound straight through the pipeline as-is."
                                activeColor="bg-violet-950/30 border-violet-700/50"
                            />
                            <RunModeCard
                                checked={toxicityReportOnly}
                                onChange={onToxicityReportOnlyChange}
                                icon={<FlaskConical className="w-3.5 h-3.5 text-red-400" />}
                                title="Toxicity Report Only"
                                desc="Forces ADMET on all compounds and skips PDB/docking entirely. Fast profiling mode."
                                activeColor="bg-red-950/30 border-red-700/50"
                            />
                        </div>
                    </section>

                    <SectionDivider />

                    {/* SECTION 2 — Analogue Generation */}
                    <section className={analoguesOff ? "opacity-40 pointer-events-none select-none" : ""}>
                        <SectionHeader
                            icon={<Layers className="w-4 h-4" />}
                            label="Analogue Generation"
                            aside={analoguesOff && (
                                <span className="text-[11px] text-violet-400 font-medium">Disabled — Direct Score active</span>
                            )}
                        />
                        <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                            {/* Number of analogues */}
                            <div>
                                <FieldLabel hint="More analogues improves coverage but increases compute time.">
                                    Number of Analogues
                                </FieldLabel>
                                <div className="flex flex-wrap gap-1.5 mb-3">
                                    {ANALOGUE_PRESETS.map(val => (
                                        <button
                                            key={val}
                                            type="button"
                                            onClick={() => onNumAnaloguesChange(val)}
                                            className={`px-3.5 py-1.5 rounded-lg text-xs font-semibold transition-all ${numAnalogues === val
                                                ? "bg-teal-600 text-white shadow-sm shadow-teal-900/50"
                                                : "bg-gray-800 text-gray-400 hover:bg-gray-700 hover:text-gray-200"
                                                }`}
                                        >
                                            {val}
                                        </button>
                                    ))}
                                </div>
                                <input
                                    type="number"
                                    min={1}
                                    max={10000}
                                    value={numAnalogues}
                                    onChange={e => onNumAnaloguesChange(parseInt(e.target.value) || 25)}
                                    placeholder="Custom"
                                    className="input w-full text-sm"
                                />
                                <p className="text-[11px] text-gray-600 mt-2">Max 10,000 — larger batches increase compute time.</p>
                            </div>

                            {/* Scaffold lock */}
                            <div>
                                <FieldLabel hint="Every generated analogue must contain this substructure. Leave blank for no lock.">
                                    Lock Core Scaffold <span className="ml-1 text-gray-600 font-normal">SMARTS</span>
                                </FieldLabel>
                                <div className="relative mb-2">
                                    <input
                                        type="text"
                                        value={lockedScaffoldSmarts}
                                        onChange={e => onLockedScaffoldSmartsChange(e.target.value)}
                                        placeholder="e.g. c1ccc2cc1C(=O)c1ccccc1N2"
                                        className="input font-mono text-xs w-full pr-8"
                                    />
                                    {lockedScaffoldSmarts && (
                                        <button
                                            type="button"
                                            onClick={() => onLockedScaffoldSmartsChange("")}
                                            className="absolute right-2.5 top-1/2 -translate-y-1/2 text-gray-500 hover:text-red-400 transition-colors"
                                        >
                                            <X className="w-3.5 h-3.5" />
                                        </button>
                                    )}
                                </div>
                                <div className="relative">
                                    <button
                                        type="button"
                                        onClick={() => setScaffoldPresetOpen(!scaffoldPresetOpen)}
                                        className="flex items-center gap-1.5 text-xs text-teal-500 hover:text-teal-400 transition-colors font-medium"
                                    >
                                        Common presets...
                                        <ChevronDown className={`w-3 h-3 transition-transform ${scaffoldPresetOpen ? "rotate-180" : ""}`} />
                                    </button>
                                    {scaffoldPresetOpen && (
                                        <div className="absolute top-full left-0 mt-2 w-full max-w-sm bg-gray-900 border border-gray-700/80 rounded-2xl shadow-2xl shadow-black/40 z-20 overflow-hidden">
                                            <div className="p-1.5 space-y-0.5">
                                                {SCAFFOLD_PRESETS.map(p => (
                                                    <button
                                                        key={p.label}
                                                        type="button"
                                                        onClick={() => { onLockedScaffoldSmartsChange(p.smarts); setScaffoldPresetOpen(false); }}
                                                        className="w-full text-left px-3.5 py-3 rounded-xl hover:bg-gray-800/70 transition-colors group"
                                                    >
                                                        <p className="text-sm font-semibold text-gray-200">{p.label}</p>
                                                        <p className="text-xs text-gray-500 mt-0.5">{p.hint}</p>
                                                        <p className="text-[11px] font-mono text-gray-700 group-hover:text-teal-500 transition-colors mt-1 truncate">{p.smarts}</p>
                                                    </button>
                                                ))}
                                            </div>
                                        </div>
                                    )}
                                </div>
                            </div>
                        </div>
                    </section>

                    <SectionDivider />

                    {/* SECTION 3 — Drug-likeness */}
                    <section>
                        <SectionHeader icon={<Activity className="w-4 h-4" />} label="Drug-likeness Limits" />
                        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                            {/* MW range */}
                            <div>
                                <FieldLabel hint="Lipinski RO5 default is 200–500 Da.">Molecular Weight Range (Da)</FieldLabel>
                                <div className="flex items-center gap-2">
                                    <input
                                        type="number"
                                        value={mwMin}
                                        onChange={e => onMwMinChange(Number(e.target.value))}
                                        className="input w-full text-sm"
                                        placeholder="Min"
                                    />
                                    <span className="text-gray-600 shrink-0 text-sm">–</span>
                                    <input
                                        type="number"
                                        value={mwMax}
                                        onChange={e => onMwMaxChange(Number(e.target.value))}
                                        className="input w-full text-sm"
                                        placeholder="Max"
                                    />
                                </div>
                                <p className="text-[11px] text-gray-600 mt-2">RO5 default 200–500 Da</p>
                            </div>

                            {/* Lipinski violations */}
                            <div>
                                <FieldLabel hint="Number of Lipinski RO5 criteria a compound may break and still pass.">
                                    Lipinski Violations
                                </FieldLabel>
                                <div className="space-y-1.5">
                                    {VIOLATION_OPTIONS.map(opt => (
                                        <label
                                            key={opt.value}
                                            className={`flex items-center gap-3 px-3.5 py-2.5 rounded-xl border cursor-pointer transition-all ${violationVal === opt.value
                                                ? "border-teal-600/60 bg-teal-950/25"
                                                : "border-gray-800 hover:border-gray-700 hover:bg-gray-800/30"
                                                }`}
                                        >
                                            <span className={`w-4 h-4 rounded-full border-2 flex items-center justify-center shrink-0 transition-all ${violationVal === opt.value ? "border-teal-500 bg-teal-500" : "border-gray-600"}`}>
                                                {violationVal === opt.value && <span className="w-1.5 h-1.5 rounded-full bg-white" />}
                                            </span>
                                            <input type="radio" name="lipinskiViolation" value={opt.value} checked={violationVal === opt.value} onChange={() => handleViolation(opt.value)} className="sr-only" />
                                            <div>
                                                <p className={`text-xs font-semibold ${violationVal === opt.value ? "text-gray-100" : "text-gray-400"}`}>{opt.label}</p>
                                                <p className="text-[11px] text-gray-600 mt-0.5">{opt.sub}</p>
                                            </div>
                                        </label>
                                    ))}
                                </div>
                            </div>

                            {/* Solubility */}
                            <div>
                                <FieldLabel hint="logS-based ESOL filter applied at analogue generation time.">
                                    Solubility Filter (logS)
                                </FieldLabel>
                                <div className="space-y-1.5">
                                    {SOLUBILITY_OPTIONS.map(opt => (
                                        <label
                                            key={opt.value}
                                            className={`flex items-center gap-3 px-3.5 py-2.5 rounded-xl border cursor-pointer transition-all ${solubilityFilter === opt.value
                                                ? "border-teal-600/60 bg-teal-950/25"
                                                : "border-gray-800 hover:border-gray-700 hover:bg-gray-800/30"
                                                }`}
                                        >
                                            <span className={`w-4 h-4 rounded-full border-2 flex items-center justify-center shrink-0 transition-all ${solubilityFilter === opt.value ? "border-teal-500 bg-teal-500" : "border-gray-600"}`}>
                                                {solubilityFilter === opt.value && <span className="w-1.5 h-1.5 rounded-full bg-white" />}
                                            </span>
                                            <input type="radio" name="solubilityFilter" value={opt.value} checked={solubilityFilter === opt.value} onChange={() => onSolubilityFilterChange(opt.value)} className="sr-only" />
                                            <div>
                                                <p className={`text-xs font-semibold ${solubilityFilter === opt.value ? "text-gray-100" : "text-gray-400"}`}>{opt.label}</p>
                                                <p className="text-[11px] text-gray-600 mt-0.5">{opt.sub}</p>
                                            </div>
                                        </label>
                                    ))}
                                </div>
                            </div>
                        </div>
                    </section>

                    <SectionDivider />

                    {/* SECTION 4 — ADMET Threshold Tuning */}
                    <section>
                        <SectionHeader
                            icon={<ShieldAlert className="w-4 h-4 text-teal-500" />}
                            label="ADMET Threshold Tuning"
                        />
                        <div className="grid grid-cols-1 lg:grid-cols-[200px_1fr] gap-5">

                            {/* Preset list */}
                            <div className="space-y-2">
                                <p className="text-[11px] font-semibold uppercase tracking-widest text-gray-600 mb-3">Target Profile</p>
                                {ADMET_PRESETS.map(preset => {
                                    const isActive = admetConfig.preset === preset.value;
                                    return (
                                        <button
                                            key={preset.value}
                                            type="button"
                                            onClick={() => handleAdmetPreset(preset.value)}
                                            className={`w-full text-left flex items-center gap-3 px-3.5 py-3 rounded-2xl border transition-all ${isActive
                                                ? PRESET_ACCENT[preset.accent]
                                                : "bg-gray-900/20 border-gray-800 text-gray-400 hover:border-gray-600 hover:bg-gray-800/30 hover:text-gray-300"
                                                }`}
                                        >
                                            <span className={isActive ? "" : "text-gray-600 [&_svg]:text-gray-600"}>
                                                {preset.icon}
                                            </span>
                                            <div className="min-w-0">
                                                <p className="text-xs font-bold leading-tight">{preset.label}</p>
                                                <p className="text-[11px] mt-0.5 opacity-60 leading-snug truncate">{preset.desc}</p>
                                            </div>
                                        </button>
                                    );
                                })}
                            </div>

                            {/* Sliders panel */}
                            <div className="rounded-2xl border border-gray-700/60 bg-gray-900/30 p-5 space-y-6">

                                {/* Panel header */}
                                <div className="flex items-start justify-between gap-4">
                                    <div>
                                        <p className="text-xs text-gray-400 leading-relaxed max-w-md">
                                            Adjust decision thresholds per endpoint. <span className="text-red-400 font-medium">Hard-fail</span> removes compounds from the pipeline.{" "}
                                            <span className="text-blue-400 font-medium">Advisory</span> flags are informational only. Switch to{" "}
                                            <button
                                                type="button"
                                                onClick={() => handleAdmetPreset("custom")}
                                                className="text-amber-400 hover:text-amber-300 transition-colors font-semibold underline underline-offset-2"
                                            >
                                                Custom
                                            </button>{" "}
                                            to save manual changes.
                                        </p>
                                    </div>
                                    {admetConfig.preset !== "balanced" && (
                                        <span className={`shrink-0 text-[11px] font-bold px-2.5 py-1 rounded-lg border ${PRESET_ACCENT[ADMET_PRESETS.find(p => p.value === admetConfig.preset)?.accent ?? "teal"]}`}>
                                            {admetConfig.preset}
                                        </span>
                                    )}
                                </div>

                                {/* ── Hard-fail endpoints ── */}
                                <div>
                                    <div className="flex items-center gap-2 mb-4">
                                        <AlertTriangle className="w-3.5 h-3.5 text-red-500" />
                                        <p className="text-[11px] font-bold uppercase tracking-widest text-red-400/90">
                                            Hard-Fail Endpoints
                                        </p>
                                        <span className="text-[10px] text-gray-600">compounds removed if exceeded</span>
                                    </div>
                                    <div className="space-y-5">
                                        {hardFailEndpoints.map(cfg => {
                                            const val = cfg.defaultGetter(admetConfig);
                                            const endpointCfg = admetConfig[cfg.key as keyof ADMETTuningConfig] as { enabled?: boolean } | undefined;
                                            const isEnabled = endpointCfg?.enabled ?? true;

                                            return (
                                                <div key={cfg.key} className={`space-y-3 ${!isEnabled ? "opacity-40" : ""}`}>
                                                    <div className="flex items-start justify-between gap-4">
                                                        <div>
                                                            <div className="flex items-center gap-2">
                                                                <p className="text-sm font-bold text-gray-100">{cfg.label}</p>
                                                                <button
                                                                    type="button"
                                                                    onClick={() => handleAdmetToggle(cfg.key)}
                                                                    className={`text-[10px] font-semibold px-1.5 py-0.5 rounded border transition-colors ${isEnabled
                                                                        ? "border-red-900/60 text-red-400 hover:bg-red-950/30"
                                                                        : "border-gray-700 text-gray-600 hover:border-gray-600 hover:text-gray-400"
                                                                        }`}
                                                                >
                                                                    {isEnabled ? "active" : "off"}
                                                                </button>
                                                            </div>
                                                            <p className="text-[11px] text-gray-500 mt-0.5 leading-relaxed">
                                                                {cfg.subtitle} · <span className="text-teal-500">{cfg.safeNote}</span>
                                                            </p>
                                                        </div>
                                                        <div className="shrink-0 text-right">
                                                            <p className="text-[11px] text-gray-600 mb-1">{cfg.failLabel}</p>
                                                            <span className="text-sm font-mono font-bold text-red-400 bg-red-950/50 border border-red-900/50 px-2.5 py-1 rounded-lg">
                                                                {cfg.formatVal(val)}
                                                            </span>
                                                        </div>
                                                    </div>
                                                    <div className="space-y-1">
                                                        <input
                                                            type="range"
                                                            min={cfg.min}
                                                            max={cfg.max}
                                                            step={cfg.step}
                                                            value={val}
                                                            disabled={!isEnabled}
                                                            onChange={e => handleAdmetSlider(cfg.key, cfg.sliderField, parseFloat(e.target.value))}
                                                            style={sliderFillStyle(val, cfg.min, cfg.max, cfg.color, cfg.thumbColor)}
                                                            className="w-full h-2 rounded-full appearance-none cursor-pointer disabled:cursor-not-allowed
                                                                [&::-webkit-slider-thumb]:appearance-none
                                                                [&::-webkit-slider-thumb]:w-4 [&::-webkit-slider-thumb]:h-4
                                                                [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:bg-white
                                                                [&::-webkit-slider-thumb]:shadow-md [&::-webkit-slider-thumb]:shadow-black/40
                                                                [&::-webkit-slider-thumb]:border-2 [&::-webkit-slider-thumb]:border-[var(--thumb-border)]"
                                                        />
                                                        <div className="flex justify-between text-[11px] text-gray-600">
                                                            <span>{cfg.min.toFixed(2)} Strict</span>
                                                            <span>{cfg.max.toFixed(2)} Relaxed</span>
                                                        </div>
                                                    </div>
                                                </div>
                                            );
                                        })}
                                    </div>
                                </div>

                                <div className="h-px bg-gray-800" />

                                {/* ── Advisory endpoints ── */}
                                <div>
                                    <div className="flex items-center gap-2 mb-4">
                                        <Info className="w-3.5 h-3.5 text-blue-500" />
                                        <p className="text-[11px] font-bold uppercase tracking-widest text-blue-400/90">
                                            Advisory Endpoints
                                        </p>
                                        <span className="text-[10px] text-gray-600">flagged but not removed</span>
                                    </div>
                                    <div className="space-y-5">
                                        {advisoryEndpoints.map(cfg => {
                                            const val = cfg.defaultGetter(admetConfig);
                                            const endpointCfg = admetConfig[cfg.key as keyof ADMETTuningConfig] as { enabled?: boolean } | undefined;
                                            const isEnabled = endpointCfg?.enabled ?? true;
                                            const isBBB = cfg.key === "bbb_penetration";

                                            return (
                                                <div key={cfg.key} className={`space-y-3 ${!isEnabled ? "opacity-40" : ""}`}>
                                                    <div className="flex items-start justify-between gap-4">
                                                        <div>
                                                            <div className="flex items-center gap-2">
                                                                <p className="text-sm font-bold text-gray-100">{cfg.label}</p>
                                                                <button
                                                                    type="button"
                                                                    onClick={() => handleAdmetToggle(cfg.key)}
                                                                    className={`text-[10px] font-semibold px-1.5 py-0.5 rounded border transition-colors ${isEnabled
                                                                        ? isBBB
                                                                            ? "border-violet-900/60 text-violet-400 hover:bg-violet-950/30"
                                                                            : "border-blue-900/60 text-blue-400 hover:bg-blue-950/30"
                                                                        : "border-gray-700 text-gray-600 hover:border-gray-600 hover:text-gray-400"
                                                                        }`}
                                                                >
                                                                    {isEnabled ? "active" : "off"}
                                                                </button>
                                                            </div>
                                                            <p className="text-[11px] text-gray-500 mt-0.5 leading-relaxed">
                                                                {cfg.subtitle} · <span className={isBBB ? "text-violet-400" : "text-teal-500"}>{cfg.safeNote}</span>
                                                            </p>
                                                        </div>
                                                        <div className="shrink-0 text-right">
                                                            <p className="text-[11px] text-gray-600 mb-1">{cfg.failLabel}</p>
                                                            <span className={`text-sm font-mono font-bold px-2.5 py-1 rounded-lg border ${isBBB
                                                                ? "text-violet-300 bg-violet-950/50 border-violet-900/50"
                                                                : "text-blue-300 bg-blue-950/50 border-blue-900/50"
                                                                }`}>
                                                                {cfg.formatVal(val)}
                                                            </span>
                                                        </div>
                                                    </div>
                                                    <div className="space-y-1">
                                                        <input
                                                            type="range"
                                                            min={cfg.min}
                                                            max={cfg.max}
                                                            step={cfg.step}
                                                            value={val}
                                                            disabled={!isEnabled}
                                                            onChange={e => handleAdmetSlider(cfg.key, cfg.sliderField, parseFloat(e.target.value))}
                                                            style={sliderFillStyle(val, cfg.min, cfg.max, cfg.color, cfg.thumbColor)}
                                                            className="w-full h-2 rounded-full appearance-none cursor-pointer disabled:cursor-not-allowed
                                                                [&::-webkit-slider-thumb]:appearance-none
                                                                [&::-webkit-slider-thumb]:w-4 [&::-webkit-slider-thumb]:h-4
                                                                [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:bg-white
                                                                [&::-webkit-slider-thumb]:shadow-md [&::-webkit-slider-thumb]:shadow-black/40
                                                                [&::-webkit-slider-thumb]:border-2 [&::-webkit-slider-thumb]:border-[var(--thumb-border)]"
                                                        />
                                                        <div className="flex justify-between text-[11px] text-gray-600">
                                                            <span>{cfg.min.toFixed(2)}</span>
                                                            <span>{cfg.max.toFixed(2)}</span>
                                                        </div>
                                                    </div>
                                                </div>
                                            );
                                        })}
                                    </div>
                                </div>
                            </div>
                        </div>
                    </section>

                    <SectionDivider />

                    {/* SECTION 5 — Pipeline Steps */}
                    <section>
                        <SectionHeader
                            icon={<Layers className="w-4 h-4" />}
                            label="Active Pipeline Steps"
                            aside={toxOnly && (
                                <span className="text-[11px] text-red-400 font-medium">Overridden by Toxicity Report mode</span>
                            )}
                        />
                        <div className={`grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3 ${toxOnly ? "opacity-40 pointer-events-none select-none" : ""}`}>
                            {STEPS.map(step => {
                                let isChecked = pipelineSteps[step.key];
                                if (toxOnly) {
                                    if (step.key === "admet") isChecked = true;
                                    else if (step.key !== "drug_likeness") isChecked = false;
                                }
                                return (
                                    <label
                                        key={step.key}
                                        className={`flex items-start gap-3.5 p-4 rounded-2xl border transition-all ${isChecked
                                            ? "bg-gray-800/50 border-gray-700"
                                            : "bg-transparent border-gray-800/50 hover:border-gray-700 hover:bg-gray-800/20"
                                            } ${step.locked ? "opacity-60 cursor-not-allowed" : "cursor-pointer"}`}
                                    >
                                        <div className="pt-0.5">
                                            <span className={`flex w-4 h-4 rounded-md border-2 items-center justify-center transition-all ${isChecked ? "bg-teal-600 border-teal-600" : "bg-gray-800 border-gray-600"}`}>
                                                {isChecked && (
                                                    <svg className="w-2.5 h-2.5 text-white" viewBox="0 0 10 10" fill="none">
                                                        <path d="M2 5l2.5 2.5 3.5-4" stroke="currentColor" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" />
                                                    </svg>
                                                )}
                                            </span>
                                        </div>
                                        <input
                                            type="checkbox"
                                            disabled={step.locked || toxOnly}
                                            checked={isChecked}
                                            onChange={() => toggleStep(step.key)}
                                            className="sr-only"
                                        />
                                        <div>
                                            <p className={`text-sm font-semibold leading-tight flex items-center gap-1.5 ${isChecked ? "text-gray-100" : "text-gray-400"}`}>
                                                <span className={isChecked ? "text-teal-500" : "text-gray-600"}>{step.icon}</span>
                                                {step.label}
                                                {step.locked && (
                                                    <span className="text-[10px] text-gray-600 font-bold uppercase tracking-wide border border-gray-700 rounded px-1 py-0.5 ml-1">
                                                        Required
                                                    </span>
                                                )}
                                            </p>
                                            <p className="text-[11px] text-gray-600 mt-1.5 leading-relaxed">{step.desc}</p>
                                        </div>
                                    </label>
                                );
                            })}
                        </div>
                    </section>

                    <SectionDivider />

                    {/* SECTION 6 — Docking Settings */}
                    <section className={dockingDisabled ? "opacity-40 pointer-events-none select-none" : ""}>
                        <SectionHeader
                            icon={<Target className="w-4 h-4" />}
                            label="Docking Settings"
                            aside={dockingDisabled && (
                                <span className="text-[11px] text-gray-500 font-medium">
                                    {toxOnly ? "Disabled — Tox-only mode" : "Disabled — Docking step off"}
                                </span>
                            )}
                        />
                        <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                            {/* Compounds to dock */}
                            <div>
                                <FieldLabel hint="Number of compounds from the ML pre-filter sent to AutoDock Vina.">
                                    Max Compounds to Dock (Top N)
                                </FieldLabel>
                                <div className="flex flex-wrap gap-1.5 mb-4">
                                    {DOCKING_COUNT_PRESETS.map(val => (
                                        <button
                                            key={val}
                                            type="button"
                                            onClick={() => { setDockingCountCustom(false); onMaxDockingCompoundsChange(val); }}
                                            className={`px-3.5 py-1.5 rounded-lg text-xs font-semibold transition-all ${!dockingCountCustom && maxDockingCompounds === val
                                                ? "bg-teal-600 text-white shadow-sm shadow-teal-900/50"
                                                : "bg-gray-800 text-gray-400 hover:bg-gray-700 hover:text-gray-200"
                                                }`}
                                        >
                                            Top {val}
                                        </button>
                                    ))}
                                    <button
                                        type="button"
                                        onClick={() => setDockingCountCustom(true)}
                                        className={`px-3.5 py-1.5 rounded-lg text-xs font-semibold transition-all ${dockingCountCustom || isCustomCount
                                            ? "bg-teal-600 text-white shadow-sm shadow-teal-900/50"
                                            : "bg-gray-800 text-gray-400 hover:bg-gray-700 hover:text-gray-200"
                                            }`}
                                    >
                                        Custom
                                    </button>
                                </div>
                                {(dockingCountCustom || isCustomCount) && (
                                    <div className="flex items-center gap-4 mb-4">
                                        <input
                                            type="range"
                                            min={1}
                                            max={50}
                                            value={maxDockingCompounds}
                                            onChange={e => onMaxDockingCompoundsChange(parseInt(e.target.value))}
                                            style={sliderFillStyle(maxDockingCompounds, 1, 50, "rgb(20 184 166 / 0.7)", "#14b8a6")}
                                            className="flex-1 h-2 rounded-full appearance-none cursor-pointer
                                                [&::-webkit-slider-thumb]:appearance-none
                                                [&::-webkit-slider-thumb]:w-4 [&::-webkit-slider-thumb]:h-4
                                                [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:bg-white
                                                [&::-webkit-slider-thumb]:shadow-md [&::-webkit-slider-thumb]:shadow-black/40
                                                [&::-webkit-slider-thumb]:border-2 [&::-webkit-slider-thumb]:border-[var(--thumb-border)]"
                                        />
                                        <span className="text-sm font-mono font-bold text-teal-400 w-8 text-right shrink-0">
                                            {maxDockingCompounds}
                                        </span>
                                    </div>
                                )}
                                {/* Compute estimate */}
                                <div className={`flex items-start gap-3 px-4 py-3 rounded-xl border text-xs ${dockTimeSev === "heavy"
                                    ? "border-red-900/50 bg-red-950/20 text-red-400"
                                    : dockTimeSev === "warn"
                                        ? "border-yellow-900/50 bg-yellow-950/20 text-yellow-400"
                                        : "border-gray-700/60 bg-gray-900/30 text-gray-400"
                                    }`}>
                                    <Clock className="w-4 h-4 shrink-0 mt-0.5" />
                                    <div>
                                        <p className="font-bold">Estimated compute: {dockTimeStr}</p>
                                        <p className="mt-0.5 opacity-70">
                                            {activeSpeed.minPerCompound} min/compound × {maxDockingCompounds} compound{maxDockingCompounds !== 1 ? "s" : ""}
                                        </p>
                                    </div>
                                </div>
                            </div>

                            {/* Vina exhaustiveness */}
                            <div>
                                <FieldLabel hint="Higher exhaustiveness increases accuracy at the cost of compute time.">
                                    Vina Exhaustiveness
                                </FieldLabel>
                                <div className="space-y-2">
                                    {DOCKING_SPEEDS.map(speed => (
                                        <label
                                            key={speed.value}
                                            className={`flex items-center justify-between px-4 py-3.5 rounded-2xl border cursor-pointer transition-all ${dockingSpeed === speed.value
                                                ? "bg-teal-950/30 border-teal-600/60 text-white"
                                                : "bg-gray-800/20 border-gray-700/50 text-gray-400 hover:bg-gray-800/40 hover:border-gray-600"
                                                }`}
                                        >
                                            <div className="flex items-center gap-3">
                                                <span className={`w-4 h-4 rounded-full border-2 flex items-center justify-center shrink-0 transition-all ${dockingSpeed === speed.value ? "border-teal-500 bg-teal-500" : "border-gray-600"}`}>
                                                    {dockingSpeed === speed.value && <span className="w-1.5 h-1.5 rounded-full bg-white" />}
                                                </span>
                                                <input
                                                    type="radio"
                                                    name="dockingSpeed"
                                                    value={speed.value}
                                                    checked={dockingSpeed === speed.value}
                                                    onChange={e => onDockingSpeedChange(e.target.value as DockingSpeed)}
                                                    className="sr-only"
                                                />
                                                <div>
                                                    <p className="text-sm font-bold">{speed.label}</p>
                                                    <p className={`text-[11px] mt-0.5 font-mono ${dockingSpeed === speed.value ? "text-teal-400/70" : "text-gray-600"}`}>
                                                        exhaustiveness{speed.tag}
                                                    </p>
                                                </div>
                                            </div>
                                            <span className={`text-xs font-semibold ${dockingSpeed === speed.value ? "text-teal-300" : "text-gray-600"}`}>
                                                {speed.detail}
                                            </span>
                                        </label>
                                    ))}
                                </div>
                            </div>
                        </div>
                    </section>

                    <SectionDivider />

                    {/* SECTION 7 — Binding Site */}
                    <section className={dockingDisabled ? "opacity-40 pointer-events-none select-none" : ""}>
                        <SectionHeader
                            icon={<MapPin className="w-4 h-4" />}
                            label="Binding Site Detection"
                            aside={bindingSiteMode !== "auto" && !dockingDisabled ? (
                                <span className="text-[11px] text-sky-400 font-medium">
                                    {bindingSiteMode === "coordinates" ? "Manual Coords active" : "Residue Range active"}
                                </span>
                            ) : undefined}
                        />
                        {/* Mode selector */}
                        <div className="flex flex-wrap gap-2 mb-6">
                            {BINDING_MODES.map(mode => (
                                <button
                                    key={mode.value}
                                    type="button"
                                    onClick={() => onBindingSiteModeChange(mode.value)}
                                    title={mode.hint}
                                    className={`px-4 py-2 rounded-xl text-xs font-semibold border transition-all ${bindingSiteMode === mode.value
                                        ? "bg-sky-600 text-white border-sky-600 shadow-sm shadow-sky-900/50"
                                        : "bg-gray-800/40 text-gray-400 border-gray-700/50 hover:bg-gray-800 hover:text-gray-200 hover:border-gray-600"
                                        }`}
                                >
                                    {mode.label}
                                </button>
                            ))}
                        </div>

                        {/* Auto mode info */}
                        {bindingSiteMode === "auto" && (
                            <div className="flex items-start gap-3 px-4 py-3.5 rounded-xl border border-sky-700/30 bg-sky-950/20">
                                <Info className="w-4 h-4 shrink-0 mt-0.5 text-sky-500" />
                                <p className="text-xs text-sky-300 leading-relaxed">
                                    The backend will detect the binding site automatically from the co-crystallised ligand in the PDB file.
                                    Use <strong className="text-sky-200">Manual Coords</strong> or <strong className="text-sky-200">Residue Range</strong> for apo structures or if auto-detection fails.
                                </p>
                            </div>
                        )}

                        {/* Coordinates mode */}
                        {bindingSiteMode === "coordinates" && (
                            <div className="rounded-2xl border border-gray-700/60 bg-gray-900/30 p-5">
                                <p className="text-[11px] font-semibold uppercase tracking-widest text-gray-600 mb-4">Pocket Centre + Search Box</p>
                                <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                                    {(["x", "y", "z"] as const).map(axis => (
                                        <div key={axis}>
                                            <FieldLabel>Centre {axis.toUpperCase()}</FieldLabel>
                                            <input
                                                type="number"
                                                step={0.1}
                                                value={bindingSiteCoords[axis]}
                                                onChange={e => onBindingSiteCoordsChange({ ...bindingSiteCoords, [axis]: parseFloat(e.target.value) || 0 })}
                                                className="input w-full text-sm font-mono"
                                                placeholder="0.0"
                                            />
                                        </div>
                                    ))}
                                    <div>
                                        <FieldLabel hint="Search box side length in Angstroms. 20 covers most binding pockets — increase for large/allosteric sites.">
                                            Box Size
                                        </FieldLabel>
                                        <input
                                            type="number"
                                            step={1}
                                            min={10}
                                            max={60}
                                            value={bindingSiteCoords.box_size}
                                            onChange={e => onBindingSiteCoordsChange({ ...bindingSiteCoords, box_size: parseFloat(e.target.value) || 20 })}
                                            className="input w-full text-sm font-mono"
                                            placeholder="20"
                                        />
                                    </div>
                                </div>
                                <p className="text-[11px] text-gray-600 mt-4">
                                    Tip: Open the PDB in PyMOL or RCSB viewer, select the ligand, and read off the XYZ coords from the centre of mass.
                                </p>
                            </div>
                        )}

                        {/* Residues mode */}
                        {bindingSiteMode === "residues" && (
                            <div className="rounded-2xl border border-gray-700/60 bg-gray-900/30 p-5">
                                <p className="text-[11px] font-semibold uppercase tracking-widest text-gray-600 mb-4">Chain + Residue Range</p>
                                <div className="grid grid-cols-1 md:grid-cols-3 gap-5">
                                    <div>
                                        <FieldLabel hint="PDB chain identifier — usually A, B, C. Case-insensitive, stored as uppercase.">Chain ID</FieldLabel>
                                        <input
                                            type="text"
                                            maxLength={2}
                                            value={bindingSiteResidues.chain}
                                            onChange={e => onBindingSiteResiduesChange({ ...bindingSiteResidues, chain: e.target.value.toUpperCase() })}
                                            className="input w-full text-sm font-mono uppercase"
                                            placeholder="A"
                                        />
                                    </div>
                                    <div>
                                        <FieldLabel hint="First residue number in the binding pocket. Use PyMOL residue IDs from the ATOM records.">Residue Start</FieldLabel>
                                        <input
                                            type="number"
                                            min={1}
                                            value={bindingSiteResidues.residue_start}
                                            onChange={e => onBindingSiteResiduesChange({ ...bindingSiteResidues, residue_start: parseInt(e.target.value) || 1 })}
                                            className="input w-full text-sm font-mono"
                                            placeholder="100"
                                        />
                                    </div>
                                    <div>
                                        <FieldLabel hint="Last residue number in the binding pocket range.">Residue End</FieldLabel>
                                        <input
                                            type="number"
                                            min={1}
                                            value={bindingSiteResidues.residue_end}
                                            onChange={e => onBindingSiteResiduesChange({ ...bindingSiteResidues, residue_end: parseInt(e.target.value) || 1 })}
                                            className="input w-full text-sm font-mono"
                                            placeholder="200"
                                        />
                                    </div>
                                </div>
                                {bindingSiteResidues.residue_end > 0 && bindingSiteResidues.residue_start > 0 && bindingSiteResidues.residue_end < bindingSiteResidues.residue_start && (
                                    <p className="text-[11px] text-red-400 mt-3 flex items-center gap-1.5">
                                        <AlertTriangle className="w-3 h-3 shrink-0" />
                                        Residue End must be ≥ Residue Start.
                                    </p>
                                )}
                                <p className="text-[11px] text-gray-600 mt-4">
                                    Tip: The backend will compute the geometric centre of all Cα atoms in this range and use it as the docking box centre.
                                </p>
                            </div>
                        )}
                    </section>

                </div>
            )}
        </div>
    );
}
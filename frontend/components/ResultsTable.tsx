"use client";

import { useState } from "react";
import {
    Download, ArrowUpDown, Trophy, BookOpen, X,
    GitBranch, FlaskConical, TestTube, Thermometer,
    Clock, Percent, ArrowRight, AlertCircle, FileText,
    ChevronDown, ChevronUp, Copy, CheckCheck, Beaker,
    Layers,
} from "lucide-react";
import type { JobResultsResponse, CompoundResult, ADMETFlagDetail } from "@/lib/api";
import { getScoreColor, getAffinityColor } from "@/lib/api";
import MoleculeCard from "./MoleculeCard";


// ── Extended retrosynthesis types (same shape as MoleculeCard) ─────────────

interface ExtendedRetrosynthesisStep {
    step_number: number;
    reaction_smarts: string;
    reaction_name?: string | null;
    reaction_type?: string | null;
    starting_materials: string[];
    product_smiles?: string | null;
    schematic?: string | null;
    reagents?: string[];
    solvents?: string[];
    conditions?: string | null;
    temperature?: string | null;
    duration?: string | null;
    yield_estimate?: string | null;
    protocol_text?: string | null;
    confidence: number;
}

interface ExtendedRetrosynthesisResult {
    feasible: boolean;
    num_steps: number;
    sa_score?: number;
    route: ExtendedRetrosynthesisStep[];
    complexity_score: number;
    difficulty_label?: string;
    synthesis_summary?: string | null;
    estimated_total_yield?: string | null;
}


// ── Constants ─────────────────────────────────────────────────────────────────

const REACTION_TYPE_LABELS: Record<string, string> = {
    transition_metal_catalysis: "Transition Metal Catalysis",
    condensation: "Condensation",
    condensation_reduction: "Condensation / Reduction",
    nucleophilic_substitution: "Nucleophilic Substitution",
    cyclisation: "Cyclisation",
    multicomponent: "Multicomponent",
    heteroatom_coupling: "Heteroatom Coupling",
    carbon_coupling: "C–C Coupling",
    fgi: "Functional Group Interconversion",
};

const DIFFICULTY_COLOR: Record<string, string> = {
    easy: "text-emerald-400",
    moderate: "text-yellow-400",
    hard: "text-orange-400",
    infeasible: "text-red-400",
};

const SCORING_GUIDE = [
    {
        metric: "Overall Score (0–100)",
        description: "A weighted composite of binding strength, safety, drug-likeness, and synthesizability. Weights: Docking 50 pts · ADMET Safety 25 pts · Drug-likeness 15 pts · Synthesis Ease 10 pts. When a pipeline step is disabled, its weight is redistributed proportionally to active steps.",
        ranges: [
            { range: "80–100", label: "Exceptional", color: "text-emerald-400", note: "Ready for in-vitro testing" },
            { range: "60–79", label: "Strong", color: "text-yellow-400", note: "Minor optimization needed" },
            { range: "40–59", label: "Moderate", color: "text-orange-400", note: "Significant issues present" },
            { range: "0–39", label: "Poor", color: "text-red-400", note: "Fails critical drug-like criteria" },
        ],
    },
    {
        metric: "Docking Affinity (kcal/mol) — 50 pts",
        description: "How tightly the drug molecule binds to the target protein's active site. More negative = stronger binding. Calculated using AutoDock Vina molecular simulation. Normalised from −3 (no binding) to −12 kcal/mol (exceptional).",
        ranges: [
            { range: "≤ −9.0", label: "Outstanding", color: "text-emerald-400", note: "Exceptionally strong binder" },
            { range: "−7.0 to −8.9", label: "Strong", color: "text-emerald-400", note: "Good clinical candidate" },
            { range: "−5.0 to −6.9", label: "Moderate", color: "text-yellow-400", note: "Needs optimization" },
            { range: "> −5.0", label: "Weak", color: "text-red-400", note: "Poor to no binding" },
        ],
    },
    {
        metric: "ADMET Safety — 25 pts",
        description: "Absorption, Distribution, Metabolism, Excretion, Toxicity profile. Starts from a risk-based base score (low=100%, moderate=55%, high=0%) then deducts per-flag penalties. hERG and Hepatotoxicity flags apply a ×2 multiplier due to clinical severity.",
        ranges: [
            { range: "Low risk, 0 flags", label: "Full score", color: "text-emerald-400", note: "25 pts" },
            { range: "Moderate risk", label: "Reduced", color: "text-yellow-400", note: "~14 pts base" },
            { range: "High risk", label: "Zero", color: "text-red-400", note: "0 pts — automatic fail" },
        ],
    },
    {
        metric: "Drug-likeness — 15 pts",
        description: "Weighted combination: Lipinski violations (50%) + LogS aqueous solubility (30%) + LogP optimality (20%). A molecule needs to dissolve and cross membranes to reach its target.",
        ranges: [
            { range: "0 violations, LogS ≥ −3", label: "Excellent", color: "text-emerald-400", note: "~15 pts" },
            { range: "1 violation, LogS −3 to −4", label: "Good", color: "text-yellow-400", note: "~11 pts" },
            { range: "2+ violations, LogS < −5", label: "Poor", color: "text-red-400", note: "< 8 pts" },
        ],
    },
    {
        metric: "Synthesis Ease — 10 pts",
        description: "How feasible it is to manufacture this molecule. Combines SA Score complexity (60%) + step count penalty (40%). Infeasible routes score 0.",
        ranges: [
            { range: "SA ≤ 3.5, ≤ 2 steps", label: "Easy", color: "text-emerald-400", note: "~10 pts" },
            { range: "SA 3.5–5.0, 3–5 steps", label: "Moderate", color: "text-yellow-400", note: "~6–8 pts" },
            { range: "SA > 6.0 or infeasible", label: "Hard/Fail", color: "text-red-400", note: "0–3 pts" },
        ],
    },
    {
        metric: "Docking Affinity (kcal/mol)",
        description: "How tightly the drug molecule binds to the target protein's active site. More negative = stronger binding. Calculated using AutoDock Vina molecular simulation.",
        ranges: [
            { range: "≤ −9.0", label: "Outstanding", color: "text-emerald-400", note: "Exceptionally strong binder" },
            { range: "−7.0 to −8.9", label: "Strong", color: "text-emerald-400", note: "Good clinical candidate" },
            { range: "−5.0 to −6.9", label: "Moderate", color: "text-yellow-400", note: "Needs optimization" },
            { range: "> −5.0", label: "Weak", color: "text-red-400", note: "Poor to no binding" },
        ],
    },
    {
        metric: "Synthesis Complexity (SA Score)",
        description: "How difficult and expensive it will be to manufacture this molecule in a lab. Based on the SA Score algorithm which analyses ring systems, stereocenters, and functional groups.",
        ranges: [
            { range: "< 15", label: "Very Easy", color: "text-emerald-400", note: "Cheap, 1–2 steps, standard reagents" },
            { range: "15–25", label: "Moderate", color: "text-yellow-400", note: "Standard for most oral drugs" },
            { range: "25–40", label: "Difficult", color: "text-orange-400", note: "Many steps, expensive reagents" },
            { range: "> 40", label: "Infeasible", color: "text-red-400", note: "Practically unsynthesizable" },
        ],
    },
    {
        metric: "Aqueous Solubility (LogS)",
        description: "How well the drug dissolves in water/blood. Poor solubility is a top reason drugs fail in clinical trials. Estimated using the ESOL (Delaney) model.",
        ranges: [
            { range: "> −1.0", label: "Very High", color: "text-emerald-400", note: "Dissolves freely — perfect score" },
            { range: "−1.0 to −3.0", label: "High", color: "text-emerald-400", note: "Excellent for oral delivery" },
            { range: "−3.0 to −4.0", label: "Good", color: "text-yellow-400", note: "Adequate for most formulations" },
            { range: "−4.0 to −5.0", label: "Borderline", color: "text-orange-400", note: "May need special formulation" },
            { range: "< −5.0", label: "Poor / Insoluble", color: "text-red-400", note: "High absorption failure risk" },
        ],
    },
    {
        metric: "hERG Inhibition",
        description: "Probability the drug blocks the hERG potassium channel in the heart. This is the most common cause of fatal drug withdrawals from markets worldwide.",
        ranges: [
            { range: "< 30%", label: "Safe", color: "text-emerald-400", note: "No cardiac risk signal" },
            { range: "30–50%", label: "Caution", color: "text-yellow-400", note: "Monitor in trials" },
            { range: "> 50%", label: "High Risk", color: "text-red-400", note: "Likely cardiotoxic — avoid" },
        ],
    },
    {
        metric: "Hepatotoxicity",
        description: "Probability the drug damages liver cells (hepatocytes). Liver toxicity is the #1 reason drugs are rejected by the FDA after Phase II trials.",
        ranges: [
            { range: "< 30%", label: "Safe", color: "text-emerald-400", note: "Low liver damage risk" },
            { range: "30–50%", label: "Moderate Risk", color: "text-yellow-400", note: "Needs monitoring" },
            { range: "> 50%", label: "Toxic", color: "text-red-400", note: "Likely hepatotoxic" },
        ],
    },
    {
        metric: "Oral Bioavailability",
        description: "What percentage of the drug actually reaches the bloodstream after swallowing a pill.",
        ranges: [
            { range: "> 70%", label: "Excellent", color: "text-emerald-400", note: "Most drug reaches target" },
            { range: "30–70%", label: "Moderate", color: "text-yellow-400", note: "Acceptable for most drugs" },
            { range: "< 30%", label: "Poor", color: "text-red-400", note: "Drug barely absorbed orally" },
        ],
    },
    {
        metric: "Caco-2 Permeability",
        description: "Measures how well the drug can pass through intestinal wall cells. A key predictor of oral absorption.",
        ranges: [
            { range: "> −5.15", label: "High Permeability", color: "text-emerald-400", note: "Crosses gut wall easily" },
            { range: "< −5.15", label: "Low Permeability", color: "text-red-400", note: "Poorly absorbed from gut" },
        ],
    },
];


// ── CSV helpers ───────────────────────────────────────────────────────────────

function cell(v: string | number | boolean | null | undefined, dp?: number): string {
    if (v == null) return "";
    if (typeof v === "boolean") return v ? "true" : "false";
    if (typeof v === "number") return dp != null ? v.toFixed(dp) : String(v);
    const s = String(v);
    return /[",|\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

function flagSeverity(flags: ADMETFlagDetail[], propName: string): string {
    return flags.find((f) => f.property_name === propName)?.severity ?? "none";
}

function buildCSV(results: JobResultsResponse): string {
    const headers = [
        "rank", "smiles", "canonical_smiles", "final_score",
        "lipinski_passed", "mw_da", "logp", "hbd", "hba", "logs", "solubility_class",
        "admet_passed",
        "herg_inhibition", "herg_severity",
        "caco2_permeability", "caco2_severity",
        "bbb_penetration", "bbb_severity",
        "hepatotoxicity", "hepatotox_severity",
        "oral_bioavailability", "oral_bioavail_severity",
        "admet_total_flags", "admet_high_severity_flags",
        "admet_flag_names", "admet_flag_implications", "admet_recommendations",
        "prefilter_passed", "prefilter_affinity_kcal", "prefilter_confidence",
        "docking_passed", "docking_best_affinity_kcal", "docking_cnn_score",
        "pose1_affinity", "pose1_rmsd_lb", "pose1_rmsd_ub",
        "pose2_affinity", "pose2_rmsd_lb", "pose2_rmsd_ub",
        "pose3_affinity", "pose3_rmsd_lb", "pose3_rmsd_ub",
        "retro_feasible", "retro_num_steps", "retro_sa_score", "retro_difficulty",
        "retro_complexity_score", "retro_estimated_total_yield",
        "retro_synthesis_summary",
        // per-step columns (up to 5 steps)
        ...([1, 2, 3, 4, 5].flatMap((n) => [
            `retro_step${n}_name`,
            `retro_step${n}_type`,
            `retro_step${n}_smarts`,
            `retro_step${n}_confidence`,
            `retro_step${n}_starting_materials`,
            `retro_step${n}_reagents`,
            `retro_step${n}_solvents`,
            `retro_step${n}_conditions`,
            `retro_step${n}_temperature`,
            `retro_step${n}_duration`,
            `retro_step${n}_yield_estimate`,
            `retro_step${n}_protocol`,
        ])),
        "retro_all_starting_materials",
    ];

    const rows = results.final_ranked_compounds.map((c: CompoundResult) => {
        const lip = c.lipinski;
        const adm = c.admet;
        const pre = c.binding_prefilter;
        const dock = c.docking;
        const ret = c.retrosynthesis as ExtendedRetrosynthesisResult | null | undefined;

        const flags = adm?.flags ?? [];
        const highFlags = flags.filter((f) => f.severity === "high");
        const poses = dock?.poses ?? [];
        const steps = ret?.route ?? [];

        const poseCol = (i: number, key: "affinity_kcal" | "rmsd_lb" | "rmsd_ub") =>
            cell(poses[i]?.[key], 3);

        const stepCols = [1, 2, 3, 4, 5].flatMap((n) => {
            const s = steps[n - 1] as ExtendedRetrosynthesisStep | undefined;
            return [
                cell(s?.reaction_name),
                cell(s?.reaction_type),
                cell(s?.reaction_smarts),
                cell(s?.confidence, 3),
                s ? `"${s.starting_materials.join(" | ")}"` : "",
                s ? `"${(s.reagents ?? []).join(" | ")}"` : "",
                s ? `"${(s.solvents ?? []).join(" | ")}"` : "",
                cell(s?.conditions),
                cell(s?.temperature),
                cell(s?.duration),
                cell(s?.yield_estimate),
                cell(s?.protocol_text),
            ];
        });

        const allStartingMats = steps.flatMap((s) => s.starting_materials).join(" | ");

        return [
            cell(c.rank), cell(c.smiles), cell(c.canonical_smiles), cell(c.final_score, 4),
            cell(lip?.passed), cell(lip?.mw, 2), cell(lip?.logp, 3),
            cell(lip?.hbd), cell(lip?.hba), cell(lip?.logs, 3), cell(lip?.solubility_class),
            cell(adm?.passed),
            cell(adm?.herg_inhibition, 4), flagSeverity(flags, "hERG Inhibition"),
            cell(adm?.caco2_permeability, 4), flagSeverity(flags, "Caco-2 Permeability"),
            cell(adm?.bbb_penetration, 4), flagSeverity(flags, "BBB Penetration"),
            cell(adm?.hepatotoxicity, 4), flagSeverity(flags, "Hepatotoxicity"),
            cell(adm?.oral_bioavailability, 4), flagSeverity(flags, "Oral Bioavailability"),
            cell(flags.length), cell(highFlags.length),
            `"${flags.map((f) => f.property_name).join(" | ")}"`,
            `"${flags.map((f) => f.implication).join(" | ")}"`,
            `"${flags.map((f) => f.recommendation).join(" | ")}"`,
            cell(pre?.passed), cell(pre?.predicted_affinity_kcal, 3), cell(pre?.confidence, 3),
            cell(dock?.passed), cell(dock?.best_affinity_kcal, 3), cell(dock?.cnn_score, 4),
            poseCol(0, "affinity_kcal"), poseCol(0, "rmsd_lb"), poseCol(0, "rmsd_ub"),
            poseCol(1, "affinity_kcal"), poseCol(1, "rmsd_lb"), poseCol(1, "rmsd_ub"),
            poseCol(2, "affinity_kcal"), poseCol(2, "rmsd_lb"), poseCol(2, "rmsd_ub"),
            cell(ret?.feasible), cell(ret?.num_steps), cell(ret?.sa_score, 3),
            cell(ret?.difficulty_label), cell(ret?.complexity_score, 3),
            cell(ret?.estimated_total_yield), cell(ret?.synthesis_summary),
            ...stepCols,
            `"${allStartingMats}"`,
        ].join(",");
    });

    return [headers.join(","), ...rows].join("\n");
}

// Exports a human-readable synthesis protocol as a .txt file for a single compound
function buildProtocolTxt(compound: CompoundResult, rank: number): string {
    const retro = compound.retrosynthesis as ExtendedRetrosynthesisResult | null | undefined;
    const smiles = compound.canonical_smiles;

    const lines: string[] = [
        "═══════════════════════════════════════════════════════════",
        `CADD PIPELINE — SYNTHESIS PROTOCOL`,
        `Compound Rank #${rank}`,
        "═══════════════════════════════════════════════════════════",
        "",
        `SMILES: ${smiles}`,
        `Score:  ${compound.final_score?.toFixed(2) ?? "—"} / 100`,
        "",
    ];

    if (!retro) {
        lines.push("No retrosynthesis data available for this compound.");
        return lines.join("\n");
    }

    lines.push(
        `Feasible:           ${retro.feasible ? "Yes" : "No"}`,
        `Steps:              ${retro.num_steps}`,
        `SA Score:           ${retro.sa_score?.toFixed(3) ?? "—"}`,
        `Difficulty:         ${retro.difficulty_label ?? "unknown"}`,
        `Complexity Score:   ${retro.complexity_score.toFixed(2)}`,
    );
    if (retro.estimated_total_yield) {
        lines.push(`Est. Total Yield:   ${retro.estimated_total_yield}`);
    }
    if (retro.synthesis_summary) {
        lines.push("", "SUMMARY", "───────", retro.synthesis_summary);
    }

    for (const step of retro.route) {
        lines.push(
            "",
            `───────────────────────────────────────────────────────────`,
            `STEP ${step.step_number}${step.reaction_name ? ` — ${step.reaction_name}` : ""}`,
            `───────────────────────────────────────────────────────────`,
        );
        if (step.reaction_type) {
            lines.push(`Type:               ${REACTION_TYPE_LABELS[step.reaction_type] ?? step.reaction_type}`);
        }
        lines.push(`Confidence:         ${(step.confidence * 100).toFixed(0)}%`);
        if (step.temperature) lines.push(`Temperature:        ${step.temperature}`);
        if (step.duration) lines.push(`Duration:           ${step.duration}`);
        if (step.yield_estimate) lines.push(`Estimated Yield:    ${step.yield_estimate}`);
        if (step.conditions) lines.push("", `Conditions: ${step.conditions}`);

        if (step.starting_materials.length > 0) {
            lines.push("", "Starting Materials:");
            step.starting_materials.forEach((sm, i) => lines.push(`  ${i + 1}. ${sm}`));
        }
        if ((step.reagents?.length ?? 0) > 0) {
            lines.push("", "Reagents:");
            step.reagents!.forEach((r, i) => lines.push(`  ${i + 1}. ${r}`));
        }
        if ((step.solvents?.length ?? 0) > 0) {
            lines.push("", "Solvents:");
            step.solvents!.forEach((s, i) => lines.push(`  ${i + 1}. ${s}`));
        }
        lines.push("", `SMARTS: ${step.reaction_smarts}`);
        if (step.schematic) lines.push(`Route:  ${step.schematic}`);

        if (step.protocol_text) {
            lines.push("", "─── Lab Protocol ───────────────────────────────────────", step.protocol_text);
        }
    }

    lines.push("", "═══════════════════════════════════════════════════════════");
    lines.push("Generated by CADD Pipeline · ADMET-AI · AutoDock Vina");
    return lines.join("\n");
}


// ── Synthesis Modal ───────────────────────────────────────────────────────────

interface SynthesisModalProps {
    compound: CompoundResult;
    rank: number;
    onClose: () => void;
}

const SynthesisStepCard = ({ step }: { step: ExtendedRetrosynthesisStep }) => {
    const [open, setOpen] = useState(true);
    const [protocolOpen, setProtocolOpen] = useState(false);

    return (
        <div className="rounded-xl border border-gray-700 bg-gray-900/60 overflow-hidden">

            {/* Step header — always visible */}
            <div
                className="flex items-center justify-between gap-3 px-4 py-3 cursor-pointer hover:bg-gray-800/40 transition-colors select-none"
                onClick={() => setOpen(!open)}
            >
                <div className="flex items-center gap-3 min-w-0">
                    <span className="flex-shrink-0 w-6 h-6 rounded-full bg-teal-900/60 border border-teal-700/60 text-teal-300 text-xs font-bold flex items-center justify-center">
                        {step.step_number}
                    </span>
                    <div className="min-w-0">
                        <p className="text-sm font-medium text-gray-100 truncate">
                            {step.reaction_name ?? "Unnamed transformation"}
                        </p>
                        {step.reaction_type && (
                            <p className="text-xs text-gray-500 mt-0.5">
                                {REACTION_TYPE_LABELS[step.reaction_type] ?? step.reaction_type}
                            </p>
                        )}
                    </div>
                </div>
                <div className="flex items-center gap-3 flex-shrink-0">
                    {step.yield_estimate && (
                        <span className="flex items-center gap-1 text-xs text-emerald-400 bg-emerald-950/40 border border-emerald-900/50 rounded-full px-2 py-0.5">
                            <Percent className="w-3 h-3" />
                            {step.yield_estimate}
                        </span>
                    )}
                    <span className="text-xs text-gray-500 font-mono">
                        {(step.confidence * 100).toFixed(0)}% conf.
                    </span>
                    {open
                        ? <ChevronUp className="w-4 h-4 text-gray-600" />
                        : <ChevronDown className="w-4 h-4 text-gray-600" />
                    }
                </div>
            </div>

            {open && (
                <div className="px-4 pb-4 pt-1 space-y-4 border-t border-gray-800/80">

                    {/* Quick stats row */}
                    {(step.temperature || step.duration || step.yield_estimate) && (
                        <div className="flex flex-wrap gap-4 pt-1">
                            {step.temperature && (
                                <div className="flex items-center gap-1.5 text-xs">
                                    <Thermometer className="w-3.5 h-3.5 text-gray-600" />
                                    <span className="text-gray-500">Temp:</span>
                                    <span className="text-gray-200">{step.temperature}</span>
                                </div>
                            )}
                            {step.duration && (
                                <div className="flex items-center gap-1.5 text-xs">
                                    <Clock className="w-3.5 h-3.5 text-gray-600" />
                                    <span className="text-gray-500">Duration:</span>
                                    <span className="text-gray-200">{step.duration}</span>
                                </div>
                            )}
                            {step.yield_estimate && (
                                <div className="flex items-center gap-1.5 text-xs">
                                    <Percent className="w-3.5 h-3.5 text-gray-600" />
                                    <span className="text-gray-500">Est. yield:</span>
                                    <span className="text-emerald-400 font-medium">{step.yield_estimate}</span>
                                </div>
                            )}
                        </div>
                    )}

                    {/* Conditions */}
                    {step.conditions && (
                        <div>
                            <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-1.5 font-semibold">
                                Conditions
                            </p>
                            <p className="text-sm text-gray-400 leading-relaxed bg-gray-800/60 rounded-lg px-3 py-2 border border-gray-700/50">
                                {step.conditions}
                            </p>
                        </div>
                    )}

                    {/* Starting materials */}
                    {step.starting_materials.length > 0 && (
                        <div>
                            <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-2 font-semibold flex items-center gap-1.5">
                                <Layers className="w-3 h-3" />
                                Starting Materials
                            </p>
                            <div className="space-y-1.5">
                                {step.starting_materials.map((sm, i) => (
                                    <div key={i} className="flex items-center gap-2">
                                        <span className="text-[10px] text-gray-600 w-4 text-right flex-shrink-0">{i + 1}.</span>
                                        <code className="flex-1 font-mono text-xs text-teal-400 bg-gray-800 px-2.5 py-1.5 rounded-lg border border-gray-700 truncate" title={sm}>
                                            {sm}
                                        </code>
                                    </div>
                                ))}
                            </div>
                        </div>
                    )}

                    {/* Reagents */}
                    {(step.reagents?.length ?? 0) > 0 && (
                        <div>
                            <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-2 font-semibold flex items-center gap-1.5">
                                <FlaskConical className="w-3 h-3" />
                                Reagents
                            </p>
                            <div className="flex flex-wrap gap-1.5">
                                {step.reagents!.map((r, i) => (
                                    <span key={i} className="text-xs text-gray-300 bg-gray-800 border border-gray-700 rounded-lg px-2.5 py-1 font-mono">
                                        {r}
                                    </span>
                                ))}
                            </div>
                        </div>
                    )}

                    {/* Solvents */}
                    {(step.solvents?.length ?? 0) > 0 && (
                        <div>
                            <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-2 font-semibold flex items-center gap-1.5">
                                <TestTube className="w-3 h-3" />
                                Solvents
                            </p>
                            <div className="flex flex-wrap gap-1.5">
                                {step.solvents!.map((s, i) => (
                                    <span key={i} className="text-xs text-blue-300/80 bg-blue-950/30 border border-blue-900/50 rounded-lg px-2.5 py-1 font-mono">
                                        {s}
                                    </span>
                                ))}
                            </div>
                        </div>
                    )}

                    {/* Schematic arrow diagram */}
                    {step.schematic && (
                        <div>
                            <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-2 font-semibold">
                                Route Schematic
                            </p>
                            <div className="flex items-center flex-wrap gap-1.5 bg-gray-800/40 rounded-lg px-3 py-2 border border-gray-700/50 overflow-x-auto">
                                {step.schematic.split("→").map((part, i, arr) => (
                                    <div key={i} className="flex items-center gap-1.5 flex-shrink-0">
                                        <code className="font-mono text-xs text-teal-400 bg-gray-900 px-2 py-1 rounded border border-gray-700 max-w-[160px] truncate block" title={part.trim()}>
                                            {part.trim()}
                                        </code>
                                        {i < arr.length - 1 && (
                                            <ArrowRight className="w-3.5 h-3.5 text-gray-600 flex-shrink-0" />
                                        )}
                                    </div>
                                ))}
                            </div>
                        </div>
                    )}

                    {/* SMARTS */}
                    <div>
                        <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-1.5 font-semibold">
                            Reaction SMARTS
                        </p>
                        <code className="block font-mono text-xs text-gray-400 bg-gray-800/60 px-3 py-2 rounded-lg border border-gray-700 break-all leading-relaxed">
                            {step.reaction_smarts}
                        </code>
                    </div>

                    {/* Lab protocol — collapsible */}
                    {step.protocol_text && (
                        <div className="border border-gray-700/80 rounded-xl overflow-hidden">
                            <button
                                onClick={() => setProtocolOpen(!protocolOpen)}
                                className="w-full flex items-center justify-between gap-2 px-4 py-2.5 text-xs text-gray-400 hover:text-gray-200 hover:bg-gray-800/60 transition-colors"
                            >
                                <div className="flex items-center gap-2">
                                    <FileText className="w-3.5 h-3.5 text-gray-500" />
                                    <span className="font-semibold uppercase tracking-wider text-[11px]">
                                        Lab Protocol
                                    </span>
                                    <span className="text-[10px] text-gray-600">
                                        Step-by-step procedure
                                    </span>
                                </div>
                                {protocolOpen
                                    ? <ChevronUp className="w-3.5 h-3.5 text-gray-600" />
                                    : <ChevronDown className="w-3.5 h-3.5 text-gray-600" />
                                }
                            </button>
                            {protocolOpen && (
                                <div className="px-4 pb-4 pt-2 border-t border-gray-700/60">
                                    <pre className="text-xs text-gray-400 leading-relaxed whitespace-pre-wrap font-mono bg-gray-900/60 rounded-lg px-3 py-2.5 border border-gray-800">
                                        {step.protocol_text}
                                    </pre>
                                </div>
                            )}
                        </div>
                    )}
                </div>
            )}
        </div>
    );
};

const SynthesisModal = ({ compound, rank, onClose }: SynthesisModalProps) => {
    const [smilescopied, setSmilesCopied] = useState(false);
    const retro = compound.retrosynthesis as ExtendedRetrosynthesisResult | null | undefined;

    const diffColor = DIFFICULTY_COLOR[retro?.difficulty_label ?? "unknown"] ?? "text-gray-500";

    const copySMILES = () => {
        navigator.clipboard.writeText(compound.canonical_smiles);
        setSmilesCopied(true);
        setTimeout(() => setSmilesCopied(false), 2000);
    };

    const downloadProtocol = () => {
        const txt = buildProtocolTxt(compound, rank);
        const blob = new Blob([txt], { type: "text/plain;charset=utf-8;" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = `synthesis_rank${rank}_${compound.canonical_smiles.slice(0, 12).replace(/[^a-zA-Z0-9]/g, "")}.txt`;
        a.click();
        URL.revokeObjectURL(url);
    };

    return (
        <div
            className="fixed inset-0 z-50 flex items-start justify-center bg-black/80 backdrop-blur-sm p-4 overflow-y-auto"
            onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
        >
            <div className="relative w-full max-w-2xl bg-gray-900 border border-gray-700 rounded-2xl shadow-2xl my-8">

                {/* Modal header */}
                <div className="flex items-center justify-between px-5 py-4 border-b border-gray-800 sticky top-0 bg-gray-900 rounded-t-2xl z-10">
                    <div className="flex items-center gap-2.5 min-w-0">
                        <GitBranch className="w-4 h-4 text-teal-400 flex-shrink-0" />
                        <div className="min-w-0">
                            <p className="text-sm font-semibold text-gray-100">
                                Synthesis Route — Rank #{rank}
                            </p>
                            <p className="text-xs text-gray-500 mt-0.5">
                                {retro
                                    ? retro.feasible
                                        ? `${retro.num_steps} step${retro.num_steps !== 1 ? "s" : ""} · SA ${retro.sa_score?.toFixed(2) ?? "—"}`
                                        : "No feasible route"
                                    : "No retrosynthesis data"
                                }
                            </p>
                        </div>
                    </div>
                    <div className="flex items-center gap-2 flex-shrink-0 ml-3">
                        <button
                            onClick={downloadProtocol}
                            title="Download full synthesis protocol as .txt"
                            className="flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg text-xs font-medium border border-gray-700 bg-gray-800/60 text-gray-400 hover:bg-gray-700 hover:text-gray-200 transition-colors"
                        >
                            <Download className="w-3.5 h-3.5" />
                            Protocol .txt
                        </button>
                        <button
                            onClick={onClose}
                            className="text-gray-600 hover:text-gray-400 transition-colors p-1"
                        >
                            <X className="w-4 h-4" />
                        </button>
                    </div>
                </div>

                <div className="px-5 py-4 space-y-5">

                    {/* SMILES */}
                    <div className="flex items-center gap-2">
                        <code className="flex-1 font-mono text-xs text-emerald-400 bg-gray-800/70 px-3 py-2 rounded-lg border border-gray-700 truncate" title={compound.canonical_smiles}>
                            {compound.canonical_smiles}
                        </code>
                        <button
                            onClick={copySMILES}
                            className="flex-shrink-0 p-2 text-gray-600 hover:text-gray-400 transition-colors"
                            title="Copy SMILES"
                        >
                            {smilescopied
                                ? <CheckCheck className="w-3.5 h-3.5 text-emerald-500" />
                                : <Copy className="w-3.5 h-3.5" />
                            }
                        </button>
                    </div>

                    {/* No data fallback */}
                    {!retro && (
                        <div className="rounded-xl bg-gray-800/40 border border-gray-700 px-4 py-6 text-center">
                            <GitBranch className="w-8 h-8 text-gray-700 mx-auto mb-2" />
                            <p className="text-sm text-gray-500">No retrosynthesis data for this compound.</p>
                            <p className="text-xs text-gray-600 mt-1">
                                Enable the retrosynthesis step and re-run the pipeline.
                            </p>
                        </div>
                    )}

                    {/* Infeasible */}
                    {retro && !retro.feasible && (
                        <div className="rounded-xl bg-red-950/20 border border-red-900/50 px-4 py-4">
                            <div className="flex items-start gap-3">
                                <AlertCircle className="w-4 h-4 text-red-500 flex-shrink-0 mt-0.5" />
                                <div>
                                    <p className="text-sm font-medium text-red-400 mb-1">No feasible synthesis route</p>
                                    {retro.synthesis_summary
                                        ? <p className="text-xs text-red-300/70 leading-relaxed">{retro.synthesis_summary}</p>
                                        : <p className="text-xs text-red-300/70">SA Score too high — molecular complexity exceeds practical synthesis limits.</p>
                                    }
                                    {retro.sa_score && (
                                        <p className="text-xs text-red-500/70 mt-1.5 font-mono">
                                            SA Score: {retro.sa_score.toFixed(3)} (infeasible threshold ≈ 6.0)
                                        </p>
                                    )}
                                </div>
                            </div>
                        </div>
                    )}

                    {/* Feasible route */}
                    {retro && retro.feasible && (
                        <>
                            {/* Route stats */}
                            <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
                                {[
                                    { label: "Steps", value: retro.num_steps },
                                    { label: "SA Score", value: retro.sa_score?.toFixed(2) ?? "—", color: diffColor },
                                    { label: "Difficulty", value: retro.difficulty_label ?? "unknown", color: diffColor },
                                    { label: "Complexity", value: retro.complexity_score.toFixed(1) },
                                ].map(({ label, value, color }) => (
                                    <div key={label} className="bg-gray-800/50 rounded-xl border border-gray-700/50 px-3 py-2.5 text-center">
                                        <p className={`text-lg font-bold ${color ?? "text-gray-100"}`}>{value}</p>
                                        <p className="text-[11px] text-gray-500 mt-0.5">{label}</p>
                                    </div>
                                ))}
                            </div>

                            {retro.estimated_total_yield && (
                                <div className="flex items-center gap-2 text-sm bg-emerald-950/30 border border-emerald-900/50 rounded-xl px-4 py-2.5">
                                    <Percent className="w-4 h-4 text-emerald-500 flex-shrink-0" />
                                    <span className="text-gray-400 text-xs">Estimated total yield (multiplicative across steps):</span>
                                    <span className="text-emerald-400 font-semibold ml-auto">{retro.estimated_total_yield}</span>
                                </div>
                            )}

                            {retro.synthesis_summary && (
                                <div className="bg-gray-800/40 border border-gray-700/60 rounded-xl px-4 py-3">
                                    <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-1.5 font-semibold">
                                        Synthesis Summary
                                    </p>
                                    <p className="text-sm text-gray-400 leading-relaxed">{retro.synthesis_summary}</p>
                                </div>
                            )}

                            {/* All reagents collected view */}
                            {retro.route.some(s => (s as ExtendedRetrosynthesisStep).reagents?.length) && (
                                <div>
                                    <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-2 font-semibold flex items-center gap-1.5">
                                        <Beaker className="w-3 h-3" />
                                        All Reagents (across all steps)
                                    </p>
                                    <div className="flex flex-wrap gap-1.5">
                                        {[...new Set(
                                            retro.route.flatMap(s => (s as ExtendedRetrosynthesisStep).reagents ?? [])
                                        )].map((r, i) => (
                                            <span key={i} className="text-xs text-gray-300 bg-gray-800 border border-gray-700 rounded-lg px-2.5 py-1 font-mono">
                                                {r}
                                            </span>
                                        ))}
                                    </div>
                                </div>
                            )}

                            {/* All starting materials */}
                            {retro.route.some(s => s.starting_materials.length > 0) && (
                                <div>
                                    <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-2 font-semibold flex items-center gap-1.5">
                                        <Layers className="w-3 h-3" />
                                        All Starting Materials
                                    </p>
                                    <div className="space-y-1.5">
                                        {[...new Set(retro.route.flatMap(s => s.starting_materials))].map((sm, i) => (
                                            <code key={i} className="block font-mono text-xs text-teal-400 bg-gray-800 px-3 py-1.5 rounded-lg border border-gray-700 truncate" title={sm}>
                                                {sm}
                                            </code>
                                        ))}
                                    </div>
                                </div>
                            )}

                            {/* Step-by-step */}
                            <div>
                                <p className="text-[11px] text-gray-600 uppercase tracking-wider mb-3 font-semibold flex items-center gap-1.5">
                                    <GitBranch className="w-3 h-3" />
                                    Step-by-Step Route
                                </p>
                                <div className="space-y-3">
                                    {retro.route.map((step) => (
                                        <SynthesisStepCard
                                            key={step.step_number}
                                            step={step as ExtendedRetrosynthesisStep}
                                        />
                                    ))}
                                </div>
                            </div>
                        </>
                    )}
                </div>
            </div>
        </div>
    );
};


// ── Main component ────────────────────────────────────────────────────────────

interface ResultsTableProps {
    results: JobResultsResponse;
}

type SortKey = "rank" | "docking" | "admet" | "score";

export default function ResultsTable({ results }: ResultsTableProps) {
    const [sortKey, setSortKey] = useState<SortKey>("rank");
    const [sortAsc, setSortAsc] = useState(true);
    const [guideOpen, setGuideOpen] = useState(false);
    const [synthTarget, setSynthTarget] = useState<{ compound: CompoundResult; rank: number } | null>(null);

    const handleSort = (key: SortKey) => {
        if (sortKey === key) setSortAsc(!sortAsc);
        else { setSortKey(key); setSortAsc(true); }
    };

    const sorted = [...results.final_ranked_compounds].sort((a, b) => {
        let va = 0, vb = 0;
        if (sortKey === "rank") { va = a.rank ?? 999; vb = b.rank ?? 999; }
        if (sortKey === "score") { va = a.final_score ?? 0; vb = b.final_score ?? 0; }
        if (sortKey === "docking") { va = a.docking?.best_affinity_kcal ?? 0; vb = b.docking?.best_affinity_kcal ?? 0; }
        if (sortKey === "admet") { va = a.admet?.flags.length ?? 99; vb = b.admet?.flags.length ?? 99; }
        return sortAsc ? va - vb : vb - va;
    });

    const downloadCsv = () => {
        const csv = buildCSV(results);
        const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = `cadd_results_${results.job_id.slice(0, 8)}.csv`;
        a.click();
        URL.revokeObjectURL(url);
    };

    const SortButton = ({ label, k }: { label: string; k: SortKey }) => (
        <button
            onClick={() => handleSort(k)}
            className={`flex items-center gap-1 text-xs uppercase tracking-wider transition-colors ${sortKey === k ? "text-emerald-400" : "text-gray-500 hover:text-gray-300"
                }`}
        >
            {label}
            <ArrowUpDown className="w-3 h-3" />
        </button>
    );

    return (
        <div className="space-y-5">

            {/* Stats summary */}
            <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
                {[
                    { label: "Analogues Generated", value: results.total_analogues_generated },
                    { label: "After Drug-likeness", value: results.compounds_after_lipinski },
                    { label: "After ADMET", value: results.compounds_after_admet },
                    { label: "Docked", value: results.compounds_docked },
                ].map(({ label, value }) => (
                    <div key={label} className="card text-center">
                        <p className="text-2xl font-bold text-gray-100">{value}</p>
                        <p className="text-xs text-gray-500 mt-0.5">{label}</p>
                    </div>
                ))}
            </div>

            {/* Top compound highlight */}
            {sorted[0] && (
                <div className="p-4 rounded-xl bg-gradient-to-r from-emerald-950/40 to-teal-950/30 border border-emerald-800/40">
                    <div className="flex items-center gap-2 mb-1">
                        <Trophy className="w-4 h-4 text-yellow-400" />
                        <p className="text-sm font-medium text-gray-200">Top Candidate</p>
                        {/* ✅ FIX: score is already normalised to 100 by backend scorer */}
                        <span className={`text-sm font-bold ${getScoreColor(sorted[0].final_score ?? 0)}`}>
                            {sorted[0].final_score?.toFixed(1)} / 100
                        </span>
                    </div>
                    <div className="flex items-center gap-4 flex-wrap">
                        {sorted[0].docking && (
                            <p className="text-xs text-gray-400">
                                Best affinity:{" "}
                                <span className={`font-mono font-medium ${getAffinityColor(sorted[0].docking.best_affinity_kcal)}`}>
                                    {sorted[0].docking.best_affinity_kcal.toFixed(3)} kcal/mol
                                </span>
                            </p>
                        )}
                        {(sorted[0].retrosynthesis as ExtendedRetrosynthesisResult | null | undefined)?.feasible && (
                            <p className="text-xs text-gray-500">
                                {(sorted[0].retrosynthesis as ExtendedRetrosynthesisResult).num_steps}-step synthesis
                            </p>
                        )}
                        {(sorted[0].retrosynthesis as ExtendedRetrosynthesisResult | null | undefined)?.feasible && (
                            <button
                                onClick={() => setSynthTarget({ compound: sorted[0], rank: sorted[0].rank ?? 1 })}
                                className="flex items-center gap-1.5 text-xs text-teal-400 hover:text-teal-300 border border-teal-800/60 hover:border-teal-700 rounded-lg px-2.5 py-1 transition-colors bg-teal-950/30"
                            >
                                <GitBranch className="w-3 h-3" />
                                View synthesis route
                            </button>
                        )}
                    </div>
                </div>
            )}

            {/* Sort controls + actions */}
            <div className="flex items-center justify-between flex-wrap gap-3">
                <div className="flex items-center gap-4">
                    <span className="text-xs text-gray-600">Sort by:</span>
                    <SortButton label="Rank" k="rank" />
                    <SortButton label="Score" k="score" />
                    <SortButton label="Docking" k="docking" />
                    <SortButton label="ADMET" k="admet" />
                </div>
                <div className="flex items-center gap-2">
                    <button
                        onClick={() => setGuideOpen(true)}
                        className="btn-secondary text-xs py-1.5 px-3"
                    >
                        <BookOpen className="w-3.5 h-3.5" />
                        Scoring Guide
                    </button>
                    <button
                        onClick={downloadCsv}
                        className="btn-secondary text-xs py-1.5 px-3"
                        title="Exports all pipeline data: ADMET values + severities, docking poses, retrosynthesis routes with reagents + protocols"
                    >
                        <Download className="w-3.5 h-3.5" />
                        Export CSV
                    </button>
                </div>
            </div>

            {/* Compound cards */}
            <div className="space-y-3">
                {sorted.map((compound, i) => {
                    const retro = compound.retrosynthesis as ExtendedRetrosynthesisResult | null | undefined;
                    return (
                        <div key={compound.canonical_smiles} className="relative group">
                            <MoleculeCard
                                compound={compound}
                                jobId={results.job_id}
                                index={i}
                            />
                            {/* Synthesis button — appears on hover, top-right corner of card */}
                            {retro && (
                                <button
                                    onClick={() => setSynthTarget({ compound, rank: compound.rank ?? i + 1 })}
                                    className={`
                                        absolute top-3 right-10 z-10
                                        flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg text-xs font-medium
                                        border transition-all
                                        ${retro.feasible
                                            ? "border-teal-800/60 bg-teal-950/50 text-teal-400 hover:border-teal-600 hover:bg-teal-900/50"
                                            : "border-red-900/50 bg-red-950/30 text-red-500/70 hover:border-red-800 hover:text-red-400"
                                        }
                                        opacity-0 group-hover:opacity-100 pointer-events-none group-hover:pointer-events-auto
                                    `}
                                    title="View synthesis route and lab protocol"
                                >
                                    <GitBranch className="w-3 h-3" />
                                    {retro.feasible ? "Synthesis" : "Infeasible"}
                                </button>
                            )}
                        </div>
                    );
                })}
            </div>

            {/* Synthesis Modal */}
            {synthTarget && (
                <SynthesisModal
                    compound={synthTarget.compound}
                    rank={synthTarget.rank}
                    onClose={() => setSynthTarget(null)}
                />
            )}

            {/* Scoring Guide Modal */}
            {guideOpen && (
                <div className="fixed inset-0 z-50 flex items-start justify-center bg-black/70 backdrop-blur-sm p-4 overflow-y-auto">
                    <div className="relative w-full max-w-2xl bg-gray-900 border border-gray-700 rounded-2xl shadow-2xl my-8">
                        <div className="flex items-center justify-between p-5 border-b border-gray-800">
                            <div className="flex items-center gap-2">
                                <BookOpen className="w-4 h-4 text-emerald-400" />
                                <h2 className="text-sm font-semibold text-gray-100">Scoring Guide</h2>
                            </div>
                            <button
                                onClick={() => setGuideOpen(false)}
                                className="text-gray-600 hover:text-gray-400 transition-colors"
                            >
                                <X className="w-4 h-4" />
                            </button>
                        </div>

                        <div className="px-5 pt-4 pb-2">
                            <p className="text-xs text-gray-500 leading-relaxed">
                                This pipeline evaluates drug candidates across multiple dimensions. Weights: <span className="text-gray-400 font-medium">Docking 50 pts · ADMET Safety 25 pts · Drug-likeness 15 pts · Synthesis Ease 10 pts</span>. When a step is disabled, its weight redistributes proportionally to active steps, always summing to 100.
                            </p>
                        </div>

                        <div className="px-5 pb-5 space-y-5">
                            {SCORING_GUIDE.map((item) => (
                                <div key={item.metric} className="p-4 rounded-xl bg-gray-800/50 border border-gray-700/50">
                                    <p className="text-sm font-medium text-gray-200 mb-1">{item.metric}</p>
                                    <p className="text-xs text-gray-500 mb-3 leading-relaxed">{item.description}</p>
                                    <div className="space-y-1.5">
                                        {item.ranges.map((r) => (
                                            <div key={r.range} className="flex items-center gap-3 text-xs">
                                                <span className="font-mono text-gray-400 w-24 flex-shrink-0">{r.range}</span>
                                                <span className={`font-medium w-24 flex-shrink-0 ${r.color}`}>{r.label}</span>
                                                <span className="text-gray-600">{r.note}</span>
                                            </div>
                                        ))}
                                    </div>
                                </div>
                            ))}
                        </div>
                    </div>
                </div>
            )}
        </div>
    );
}
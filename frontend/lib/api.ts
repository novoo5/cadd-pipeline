// All communication between the Next.js frontend and the FastAPI backend lives here.
// Never call fetch() directly from a component — always go through this module.

const BACKEND_URL = "https://novoo5-cadd-backend.hf.space";

// ── Enums / union types ───────────────────────────────────────────────────────

export type DockingSpeed = "fast" | "balanced" | "thorough";
export type BindingSiteMode = "auto" | "coordinates" | "residues";
export type StepStatus = "waiting" | "running" | "done" | "skipped" | "failed";
export type JobStatus = "queued" | "running" | "done" | "failed";
export type SolubilityFilterMode = "soluble_only" | "allow_slightly" | "all";
export type DockingFileType = "pdbqt" | "pose_pdb" | "complex_pdb";
export type ADMETPreset = "balanced" | "cns" | "oral" | "custom";
export type ADMETCategory =
    | "absorption"
    | "distribution"
    | "metabolism"
    | "excretion"
    | "toxicity"
    | "physchem"
    | "other";
export type NarrativeLevel = "info" | "good" | "warning" | "danger";
export type OverallRisk = "low" | "moderate" | "high";
export type ComplexityCategory = "easy" | "moderate" | "hard" | "very hard" | "infeasible";

// ── Pipeline config types ─────────────────────────────────────────────────────

export interface PipelineSteps {
    drug_likeness: boolean;
    admet: boolean;
    binding_prefilter: boolean;
    docking: boolean;
    retrosynthesis: boolean;
}

export interface BindingSiteCoords {
    x: number;
    y: number;
    z: number;
    box_size: number;
}

export interface BindingSiteResidues {
    chain: string;
    residue_start: number;
    residue_end: number;
}

// ── ADMET config types (mirrors backend ADMETTuningConfig) ────────────────────

export interface ADMETThresholdConfig {
    enabled?: boolean;
    cutoff: number;
    direction: "above" | "below";
    severity_high: number | null;
    hard_fail: boolean;
    threshold_str?: string | null;
    implication?: string | null;
    recommendation?: string | null;
}

export interface ADMETTuningConfig {
    preset?: ADMETPreset;
    include_full_profile?: boolean;
    herg_inhibition?: ADMETThresholdConfig;
    hepatotoxicity?: ADMETThresholdConfig;
    caco2_permeability?: ADMETThresholdConfig;
    oral_bioavailability?: ADMETThresholdConfig;
    bbb_penetration?: ADMETThresholdConfig;
}

export const ADMET_PRESET_DEFAULTS: Record<ADMETPreset, Partial<ADMETTuningConfig>> = {
    balanced: {
        preset: "balanced",
        herg_inhibition: { cutoff: 0.5, direction: "above", severity_high: 0.85, hard_fail: true },
        hepatotoxicity: { cutoff: 0.5, direction: "above", severity_high: 0.85, hard_fail: true },
        caco2_permeability: { cutoff: -5.15, direction: "below", severity_high: -6.5, hard_fail: false },
        oral_bioavailability: { cutoff: 0.3, direction: "below", severity_high: 0.1, hard_fail: false },
        bbb_penetration: { cutoff: 0.3, direction: "below", severity_high: null, hard_fail: false },
    },
    cns: {
        preset: "cns",
        herg_inhibition: { cutoff: 0.5, direction: "above", severity_high: 0.85, hard_fail: true },
        hepatotoxicity: { cutoff: 0.5, direction: "above", severity_high: 0.85, hard_fail: true },
        caco2_permeability: { cutoff: -5.15, direction: "below", severity_high: -6.5, hard_fail: false },
        oral_bioavailability: { cutoff: 0.3, direction: "below", severity_high: 0.1, hard_fail: false },
        // CNS: tighter BBB requirement — compounds must penetrate
        bbb_penetration: { cutoff: 0.6, direction: "below", severity_high: 0.3, hard_fail: false },
    },
    oral: {
        preset: "oral",
        herg_inhibition: { cutoff: 0.5, direction: "above", severity_high: 0.85, hard_fail: true },
        hepatotoxicity: { cutoff: 0.5, direction: "above", severity_high: 0.85, hard_fail: true },
        // Oral: tighter permeability and bioavailability thresholds
        caco2_permeability: { cutoff: -5.0, direction: "below", severity_high: -6.0, hard_fail: false },
        oral_bioavailability: { cutoff: 0.4, direction: "below", severity_high: 0.2, hard_fail: false },
        bbb_penetration: { cutoff: 0.3, direction: "below", severity_high: null, hard_fail: false },
    },
    custom: {
        preset: "custom",
    },
};

export const ADMET_ENDPOINT_LABELS: Record<string, string> = {
    herg_inhibition: "hERG Inhibition",
    hepatotoxicity: "Hepatotoxicity",
    caco2_permeability: "Caco-2 Permeability",
    oral_bioavailability: "Oral Bioavailability",
    bbb_penetration: "BBB Penetration",
};

export const ADMET_ENDPOINT_UNITS: Record<string, string> = {
    herg_inhibition: "probability",
    hepatotoxicity: "probability",
    caco2_permeability: "log cm/s",
    oral_bioavailability: "probability",
    bbb_penetration: "probability",
};

export const ADMET_ENDPOINT_DESCRIPTIONS: Record<string, string> = {
    herg_inhibition:
        "Predicted probability of inhibiting the hERG potassium channel. High values risk QT prolongation / arrhythmia. Hard-fail if >cutoff AND >severity_high.",
    hepatotoxicity:
        "Predicted probability of drug-induced liver injury (DILI). Leading cause of FDA drug withdrawal. Hard-fail if >cutoff AND >severity_high.",
    caco2_permeability:
        "Predicted intestinal permeability (Caco-2 Wang model). More negative = lower absorption. Informational only — not all drugs require high permeability.",
    oral_bioavailability:
        "Predicted fraction of drug reaching systemic circulation orally (Bioavailability_Ma, binary classifier). Informational — ~30% FPR on diverse scaffolds.",
    bbb_penetration:
        "Predicted blood-brain barrier penetration probability (BBB Martins). Informational — only flag if CNS target is intended.",
};

// ── Job request ───────────────────────────────────────────────────────────────

export interface JobRequest {
    smiles: string;
    pdb_id?: string;
    pdb_content?: string;
    num_analogues: number;
    direct_score_only?: boolean;
    solubility_filter?: SolubilityFilterMode;
    toxicity_report_only?: boolean;
    admet_config?: ADMETTuningConfig;
    pipeline_steps: PipelineSteps;
    binding_site_mode: BindingSiteMode;
    binding_site_coords?: BindingSiteCoords;
    binding_site_residues?: BindingSiteResidues;
    docking_speed: DockingSpeed;
    max_docking_compounds?: number;
    mw_min?: number;
    mw_max?: number;
    max_lipinski_violations?: number | null;
    locked_scaffold_smarts?: string | null;
}

// ── Response types ────────────────────────────────────────────────────────────

export interface JobSubmitResponse {
    job_id: string;
    message: string;
    status_url: string;
    results_url: string;
}

export interface PipelineStepInfo {
    name: string;
    status: StepStatus;
    message?: string | null;
    duration_seconds?: number | null;
    progress_current?: number | null;
    progress_total?: number | null;
}

export interface JobStatusResponse {
    job_id: string;
    status: JobStatus;
    created_at: string;
    updated_at: string;
    steps: PipelineStepInfo[];
    error?: string;
}

// ── Pipeline result types ─────────────────────────────────────────────────────

export interface LipinskiResult {
    passed: boolean;
    mw: number;
    logp: number;
    hbd: number;
    hba: number;
    logs: number;
    solubility_class: string;
}

export interface ADMETFlagDetail {
    property_name: string;
    value: number;
    threshold: string;
    direction: "above" | "below";
    severity: "high" | "moderate" | "low";
    implication: string;
    recommendation: string;
}

export interface ADMETEndpointValue {
    key: string;
    label: string;
    category: ADMETCategory;
    value: number;
    display_value: string;
    unit?: string | null;
    interpretation?: string | null;
    threshold_applied?: string | null;
    severity?: "high" | "moderate" | "low" | null;
    hard_fail: boolean;
    triggered: boolean;
}

export interface ADMETNarrativeBlock {
    title: string;
    level: NarrativeLevel;
    body: string;
}

export interface ADMETResult {
    passed: boolean;
    overall_risk: OverallRisk;
    herg_inhibition: number;
    caco2_permeability: number;
    bbb_penetration: number;
    hepatotoxicity: number;
    oral_bioavailability: number;
    flags: ADMETFlagDetail[];
    flag_summary: string[];
    decision_basis: string[];
    narrative: ADMETNarrativeBlock[];
    endpoint_table: ADMETEndpointValue[];
    extra_properties: Record<string, number>;
    thresholds_used: Record<string, ADMETThresholdConfig>;
}

export interface BindingPrefilterResult {
    passed: boolean;
    predicted_affinity_kcal: number;
    confidence: number;
}

export interface DockingPose {
    rank: number;
    affinity_kcal: number;
    rmsd_lb: number;
    rmsd_ub: number;
}

export interface DockingResult {
    passed: boolean;
    best_affinity_kcal: number;
    cnn_score: number;
    poses: DockingPose[];
}

// ── Retrosynthesis result types (UPDATED — mirrors retrosynthesis.py) ─────────

export interface RetrosynthesisStep {
    step_number: number;

    // Identity
    reaction_name?: string | null;
    reaction_type?: string | null;
    reaction_smarts: string;

    // Materials
    starting_materials: string[];
    product_smiles?: string | null;

    // Lab protocol (rich data from retrosynthesis.py)
    schematic?: string | null;
    reagents: string[];
    conditions?: string | null;
    protocol_steps: string[];

    // Metadata
    confidence: number;
    estimated_yield?: string | null;
    notes?: string | null;
}

export interface RetrosynthesisResult {
    feasible: boolean;
    num_steps: number;
    route: RetrosynthesisStep[];
    complexity_score: number;

    // Extended fields from retrosynthesis.py
    sa_score?: number | null;
    complexity_category?: ComplexityCategory | null;
    functional_groups_detected: string[];
    disconnection_points: string[];
    synthesis_overview?: string | null;
    alternative_routes_hint?: string | null;
}

// ── Compound result ───────────────────────────────────────────────────────────

export interface CompoundResult {
    smiles: string;
    canonical_smiles: string;
    rank?: number;
    final_score?: number;
    lipinski?: LipinskiResult;
    admet?: ADMETResult;
    binding_prefilter?: BindingPrefilterResult;
    docking?: DockingResult;
    retrosynthesis?: RetrosynthesisResult;
}

export interface BindingSiteInfo {
    detection_mode: "native_ligand" | "protein_centroid" | "user_coordinates" | "user_residues";
    detected_ligand_name: string | null;
    center_x: number;
    center_y: number;
    center_z: number;
    box_size: number;
}

export interface JobResultsResponse {
    job_id: string;
    base_smiles: string;
    total_analogues_generated: number;
    compounds_after_lipinski: number;
    compounds_after_admet: number;
    compounds_after_prefilter: number;
    compounds_docked: number;
    final_ranked_compounds: CompoundResult[];
    binding_site_info?: BindingSiteInfo | null;
    created_at: string;
    completed_at: string;
}

export interface ScoreBreakdownItem {
    raw: string;
    contribution: number;
    max_possible: number;
}

export interface ScoreBreakdown {
    docking_affinity?: ScoreBreakdownItem;
    cnn_score?: ScoreBreakdownItem;
    admet_safety?: ScoreBreakdownItem;
    solubility?: ScoreBreakdownItem;
    synthesis_ease?: ScoreBreakdownItem;
    binding_prefilter?: ScoreBreakdownItem;
    mw_fragment_penalty?: boolean;
    final_score?: number;
}

// ── API functions ─────────────────────────────────────────────────────────────

export async function submitJob(request: JobRequest): Promise<JobSubmitResponse> {
    const response = await fetch(`${BACKEND_URL}/api/v1/jobs`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(request),
    });
    if (!response.ok) {
        const error = await response.json().catch(() => ({}));
        throw new Error(
            error?.detail || `Job submission failed with status ${response.status}`
        );
    }
    return response.json();
}

export async function submitJobWithFile(
    smiles: string,
    pdbFile: File | null,
    options: {
        num_analogues: number;
        docking_speed: DockingSpeed;
        max_docking_compounds?: number;
        binding_site_mode: BindingSiteMode;
        pipeline_steps: PipelineSteps;
        direct_score_only?: boolean;
        mw_min?: number;
        mw_max?: number;
        max_lipinski_violations?: number | null;
        solubility_filter?: SolubilityFilterMode;
        toxicity_report_only?: boolean;
        locked_scaffold_smarts?: string;
        admet_config?: ADMETTuningConfig;
    }
): Promise<JobSubmitResponse> {
    if (pdbFile !== null) {
        const formData = new FormData();
        formData.append("smiles", smiles);
        formData.append("pdb_file", pdbFile);
        formData.append("num_analogues", String(options.num_analogues));
        formData.append("docking_speed", options.docking_speed);
        formData.append("max_docking_compounds", String(options.max_docking_compounds ?? 10));
        formData.append("binding_site_mode", options.binding_site_mode);
        formData.append("direct_score_only", String(options.direct_score_only ?? false));
        formData.append("mw_min", String(options.mw_min ?? 200));
        formData.append("mw_max", String(options.mw_max ?? 500));
        formData.append(
            "max_lipinski_violations",
            String(
                options.max_lipinski_violations === null
                    ? -1
                    : (options.max_lipinski_violations ?? 1)
            )
        );
        formData.append("solubility_filter", options.solubility_filter ?? "all");
        formData.append(
            "toxicity_report_only",
            String(options.toxicity_report_only ?? false)
        );
        if (options.locked_scaffold_smarts) {
            formData.append("locked_scaffold_smarts", options.locked_scaffold_smarts);
        }
        if (options.admet_config) {
            formData.append("admet_config", JSON.stringify(options.admet_config));
        }

        const response = await fetch(`${BACKEND_URL}/api/v1/jobs/upload`, {
            method: "POST",
            body: formData,
        });
        if (!response.ok) {
            const error = await response.json().catch(() => ({}));
            throw new Error(
                error?.detail || `File upload failed with status ${response.status}`
            );
        }
        return response.json();
    }

    return submitJob({
        smiles,
        num_analogues: options.num_analogues,
        docking_speed: options.docking_speed,
        max_docking_compounds: options.max_docking_compounds ?? 10,
        binding_site_mode: options.binding_site_mode,
        pipeline_steps: options.pipeline_steps,
        direct_score_only: options.direct_score_only ?? false,
        mw_min: options.mw_min ?? 200,
        mw_max: options.mw_max ?? 500,
        max_lipinski_violations: options.max_lipinski_violations ?? 1,
        solubility_filter: options.solubility_filter ?? "all",
        toxicity_report_only: options.toxicity_report_only ?? false,
        locked_scaffold_smarts: options.locked_scaffold_smarts ?? null,
        admet_config: options.admet_config,
    });
}

export async function getJobStatus(jobId: string): Promise<JobStatusResponse> {
    const response = await fetch(`${BACKEND_URL}/api/v1/jobs/${jobId}/status`, {
        cache: "no-store",
    });
    if (!response.ok) {
        if (response.status === 404)
            throw new Error(`Job '${jobId}' not found. It may have expired.`);
        throw new Error(`Status fetch failed: HTTP ${response.status}`);
    }
    return response.json();
}

export async function getJobResults(jobId: string): Promise<JobResultsResponse | null> {
    const response = await fetch(`${BACKEND_URL}/api/v1/jobs/${jobId}/results`, {
        cache: "no-store",
    });
    if (response.status === 202) return null;
    if (!response.ok) {
        const error = await response.json().catch(() => ({}));
        throw new Error(error?.detail || `Results fetch failed: HTTP ${response.status}`);
    }
    return response.json();
}

export async function getScoreBreakdown(
    jobId: string,
    compoundIndex: number
): Promise<ScoreBreakdown> {
    const response = await fetch(
        `${BACKEND_URL}/api/v1/jobs/${jobId}/score-breakdown/${compoundIndex}`,
        { cache: "no-store" }
    );
    if (!response.ok)
        throw new Error(`Score breakdown fetch failed: HTTP ${response.status}`);
    return response.json();
}

export async function pingBackend(): Promise<boolean> {
    try {
        const response = await fetch(`${BACKEND_URL}/api/v1/health/ping`, {
            cache: "no-store",
            signal: AbortSignal.timeout(5000),
        });
        return response.ok;
    } catch {
        return false;
    }
}

export async function fetchPdbMetadata(pdbId: string): Promise<{
    title: string;
    resolution_angstrom: number | null;
    protein_chains: number;
} | null> {
    try {
        const query = `{ entry(entry_id: "${pdbId.toUpperCase()}") { struct { title } rcsb_entry_info { resolution_combined polymer_entity_count_protein } } }`;
        const response = await fetch("https://data.rcsb.org/graphql", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ query }),
            signal: AbortSignal.timeout(8000),
        });
        if (!response.ok) return null;
        const data = await response.json();
        const entry = data?.data?.entry;
        if (!entry) return null;
        return {
            title: entry.struct?.title || "Unknown protein",
            resolution_angstrom: entry.rcsb_entry_info?.resolution_combined?.[0] || null,
            protein_chains: entry.rcsb_entry_info?.polymer_entity_count_protein || 0,
        };
    } catch {
        return null;
    }
}

// ── Docking file downloads ────────────────────────────────────────────────────

const DOCKING_FILE_LABELS: Record<DockingFileType, string> = {
    pdbqt: "docked_poses.pdbqt",
    pose_pdb: "pose1_with_H.pdb",
    complex_pdb: "receptor_ligand_complex.pdb",
};

export function getDockingFileUrl(
    jobId: string,
    compoundIndex: number,
    type: DockingFileType
): string {
    return `${BACKEND_URL}/api/v1/jobs/${jobId}/compounds/${compoundIndex}/download/${type}`;
}

export async function downloadDockingFile(
    jobId: string,
    compoundIndex: number,
    type: DockingFileType
): Promise<void> {
    const url = getDockingFileUrl(jobId, compoundIndex, type);
    const filename = `compound_${compoundIndex}_${DOCKING_FILE_LABELS[type]}`;
    const response = await fetch(url);
    if (!response.ok) {
        const err = (await response.json().catch(() => ({}))) as { detail?: string };
        throw new Error(err.detail ?? `Download failed: HTTP ${response.status}`);
    }
    const blob = await response.blob();
    const objectUrl = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    anchor.href = objectUrl;
    anchor.download = filename;
    anchor.click();
    URL.revokeObjectURL(objectUrl);
}

// ── Helpers ───────────────────────────────────────────────────────────────────

export function getScoreColor(score: number): string {
    if (score >= 70) return "text-emerald-400";
    if (score >= 45) return "text-yellow-400";
    return "text-red-400";
}

export function getAffinityColor(affinity: number): string {
    if (affinity <= -8.0) return "text-emerald-400";
    if (affinity <= -6.0) return "text-yellow-400";
    return "text-red-400";
}

export function getFlagSeverityColor(severity: "high" | "moderate" | "low"): string {
    if (severity === "high") return "text-red-400 bg-red-950/40 border-red-800";
    if (severity === "moderate") return "text-yellow-400 bg-yellow-950/40 border-yellow-800";
    return "text-gray-400 bg-gray-800/40 border-gray-700";
}

export function getRiskBadgeColor(risk: OverallRisk): string {
    if (risk === "high") return "text-red-400 bg-red-950/40 border-red-800";
    if (risk === "moderate") return "text-yellow-400 bg-yellow-950/40 border-yellow-800";
    return "text-emerald-400 bg-emerald-950/40 border-emerald-800";
}

export function getNarrativeLevelColor(level: NarrativeLevel): string {
    if (level === "danger") return "border-red-800 bg-red-950/30 text-red-300";
    if (level === "warning") return "border-yellow-800 bg-yellow-950/30 text-yellow-300";
    if (level === "good") return "border-emerald-800 bg-emerald-950/30 text-emerald-300";
    return "border-gray-700 bg-gray-900/40 text-gray-300";
}

export function getCategoryColor(category: ADMETCategory): string {
    const map: Record<ADMETCategory, string> = {
        absorption: "text-blue-400",
        distribution: "text-purple-400",
        metabolism: "text-yellow-400",
        excretion: "text-cyan-400",
        toxicity: "text-red-400",
        physchem: "text-gray-400",
        other: "text-gray-500",
    };
    return map[category] ?? "text-gray-500";
}

export function getComplexityColor(category: ComplexityCategory | null | undefined): string {
    if (!category) return "text-gray-400";
    if (category === "easy") return "text-emerald-400";
    if (category === "moderate") return "text-yellow-400";
    if (category === "hard") return "text-orange-400";
    if (category === "very hard") return "text-red-400";
    return "text-red-600";
}

export function getSAScoreColor(score: number): string {
    if (score <= 3.5) return "text-emerald-400";
    if (score <= 5.0) return "text-yellow-400";
    if (score <= 6.5) return "text-orange-400";
    return "text-red-400";
}

export function formatProbability(value: number): string {
    return `${(value * 100).toFixed(1)}%`;
}

export function formatDuration(seconds: number): string {
    if (seconds < 60) return `${seconds.toFixed(1)}s`;
    const mins = Math.floor(seconds / 60);
    const secs = Math.round(seconds % 60);
    return `${mins}m ${secs}s`;
}
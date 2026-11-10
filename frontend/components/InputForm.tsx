"use client";

import { useState, useRef, useEffect } from "react";
import { useRouter } from "next/navigation";
import {
    Upload, FlaskConical, Loader2, AlertCircle,
    CheckCircle2, ExternalLink, X, TestTube, Zap,
    History, RotateCcw,
} from "lucide-react";
import AdvancedSettings from "./AdvancedSettings";
import {
    submitJob,
    submitJobWithFile,
    fetchPdbMetadata,
    pingBackend,
    type PipelineSteps,
    type DockingSpeed,
    type BindingSiteMode,
    type BindingSiteCoords,
    type BindingSiteResidues,
    type SolubilityFilterMode,
    type ADMETTuningConfig,
    ADMET_PRESET_DEFAULTS,
} from "@/lib/api";


const MYRICETIN_SMILES = "OC1=CC(=CC(O)=C1O)C1=C(O)C(=O)C2=C(O)C(O)=CC(O)=C2O1";

const FORM_MEMORY_KEY = "cadd_form_settings_v2"; // bumped to v2 — new admetConfig field


interface SavedFormSettings {
    pdbMode: "id" | "file";
    pdbId: string;
    numAnalogues: number;
    directScoreOnly: boolean;
    toxicityReportOnly: boolean;
    solubilityFilter: SolubilityFilterMode;
    lockedScaffoldSmarts: string;
    pipelineSteps: PipelineSteps;
    dockingSpeed: DockingSpeed;
    maxDockingCompounds: number;
    bindingSiteMode: BindingSiteMode;
    bindingSiteCoords: BindingSiteCoords;
    bindingSiteResidues: BindingSiteResidues;
    mwMin: number;
    mwMax: number;
    maxLipinskiViolations: number | null;
    admetConfig: ADMETTuningConfig; // ← NEW
}


export default function InputForm() {
    const router = useRouter();

    const [smiles, setSmiles] = useState(MYRICETIN_SMILES);
    const [pdbMode, setPdbMode] = useState<"id" | "file">("id");
    const [pdbId, setPdbId] = useState("");
    const [pdbFile, setPdbFile] = useState<File | null>(null);
    const [pdbMeta, setPdbMeta] = useState<{
        title: string;
        resolution_angstrom: number | null;
        protein_chains: number;
    } | null>(null);
    const [pdbMetaLoading, setPdbMetaLoading] = useState(false);

    const [numAnalogues, setNumAnalogues] = useState<number>(25);
    const [directScoreOnly, setDirectScoreOnly] = useState(false);
    const [toxicityReportOnly, setToxicityReportOnly] = useState(false);
    const [solubilityFilter, setSolubilityFilter] = useState<SolubilityFilterMode>("all");
    const [lockedScaffoldSmarts, setLockedScaffoldSmarts] = useState("");
    const [pipelineSteps, setPipelineSteps] = useState<PipelineSteps>({
        drug_likeness: true,
        admet: true,
        binding_prefilter: true,
        docking: true,
        retrosynthesis: true,
    });
    const [dockingSpeed, setDockingSpeed] = useState<DockingSpeed>("balanced");
    const [maxDockingCompounds, setMaxDockingCompounds] = useState<number>(10);
    const [bindingSiteMode, setBindingSiteMode] = useState<BindingSiteMode>("auto");
    const [bindingSiteCoords, setBindingSiteCoords] = useState<BindingSiteCoords>({
        x: 0, y: 0, z: 0, box_size: 20,
    });
    const [bindingSiteResidues, setBindingSiteResidues] = useState<BindingSiteResidues>({
        chain: "A", residue_start: 1, residue_end: 100,
    });
    const [mwMin, setMwMin] = useState(200);
    const [mwMax, setMwMax] = useState(500);
    const [maxLipinskiViolations, setMaxLipinskiViolations] = useState<number | null>(1);

    // ── NEW: ADMET config state ───────────────────────────────────────────────
    const [admetConfig, setAdmetConfig] = useState<ADMETTuningConfig>(
        ADMET_PRESET_DEFAULTS.balanced
    );

    const [submitting, setSubmitting] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [backendAwake, setBackendAwake] = useState<boolean | null>(null);

    const [savedAt, setSavedAt] = useState<Date | null>(null);
    const [memoryRestored, setMemoryRestored] = useState(false);

    const fileInputRef = useRef<HTMLInputElement>(null);
    const pdbDebounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);
    const saveTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);


    // ── 1. Restore ────────────────────────────────────────────────────────────
    useEffect(() => {
        try {
            const raw = localStorage.getItem(FORM_MEMORY_KEY);
            if (!raw) return;
            const s: Partial<SavedFormSettings> = JSON.parse(raw);

            if (s.pdbMode) setPdbMode(s.pdbMode);
            if (s.pdbId) setPdbId(s.pdbId);
            if (s.numAnalogues != null) setNumAnalogues(s.numAnalogues);
            if (s.directScoreOnly != null) setDirectScoreOnly(s.directScoreOnly);
            if (s.toxicityReportOnly != null) setToxicityReportOnly(s.toxicityReportOnly);
            if (s.solubilityFilter) setSolubilityFilter(s.solubilityFilter);
            if (s.lockedScaffoldSmarts != null) setLockedScaffoldSmarts(s.lockedScaffoldSmarts);
            if (s.pipelineSteps) setPipelineSteps(s.pipelineSteps);
            if (s.dockingSpeed) setDockingSpeed(s.dockingSpeed);
            if (s.maxDockingCompounds != null) setMaxDockingCompounds(s.maxDockingCompounds);
            if (s.bindingSiteMode) setBindingSiteMode(s.bindingSiteMode);
            if (s.bindingSiteCoords) setBindingSiteCoords(s.bindingSiteCoords);
            if (s.bindingSiteResidues) setBindingSiteResidues(s.bindingSiteResidues);
            if (s.mwMin != null) setMwMin(s.mwMin);
            if (s.mwMax != null) setMwMax(s.mwMax);
            if ("maxLipinskiViolations" in s) setMaxLipinskiViolations(s.maxLipinskiViolations ?? null);
            // ── NEW: restore admetConfig — fall back to balanced preset if missing
            if (s.admetConfig) setAdmetConfig(s.admetConfig);

            setMemoryRestored(true);
            setSavedAt(new Date());
        } catch { }
    }, []);


    // ── 2. Auto-save ──────────────────────────────────────────────────────────
    useEffect(() => {
        if (saveTimerRef.current) clearTimeout(saveTimerRef.current);
        saveTimerRef.current = setTimeout(() => {
            try {
                const payload: SavedFormSettings = {
                    pdbMode, pdbId, numAnalogues, directScoreOnly, toxicityReportOnly,
                    solubilityFilter, lockedScaffoldSmarts, pipelineSteps, dockingSpeed,
                    maxDockingCompounds,
                    bindingSiteMode, bindingSiteCoords, bindingSiteResidues,
                    mwMin, mwMax, maxLipinskiViolations,
                    admetConfig, // ← NEW
                };
                localStorage.setItem(FORM_MEMORY_KEY, JSON.stringify(payload));
                setSavedAt(new Date());
            } catch { }
        }, 800);
        return () => { if (saveTimerRef.current) clearTimeout(saveTimerRef.current); };
    }, [
        pdbMode, pdbId, numAnalogues, directScoreOnly, toxicityReportOnly,
        solubilityFilter, lockedScaffoldSmarts, pipelineSteps, dockingSpeed,
        maxDockingCompounds,
        bindingSiteMode, bindingSiteCoords, bindingSiteResidues,
        mwMin, mwMax, maxLipinskiViolations,
        admetConfig, // ← NEW
    ]);


    // ── 3. Clear ──────────────────────────────────────────────────────────────
    const clearSavedSettings = () => {
        try { localStorage.removeItem(FORM_MEMORY_KEY); } catch { }
        setSavedAt(null); setMemoryRestored(false);
        setPdbMode("id"); setPdbId(""); setPdbFile(null); setPdbMeta(null);
        setNumAnalogues(25); setDirectScoreOnly(false); setToxicityReportOnly(false);
        setSolubilityFilter("all"); setLockedScaffoldSmarts("");
        setPipelineSteps({
            drug_likeness: true, admet: true,
            binding_prefilter: true, docking: true, retrosynthesis: true,
        });
        setDockingSpeed("balanced");
        setMaxDockingCompounds(10);
        setBindingSiteMode("auto");
        setBindingSiteCoords({ x: 0, y: 0, z: 0, box_size: 20 });
        setBindingSiteResidues({ chain: "A", residue_start: 1, residue_end: 100 });
        setMwMin(200); setMwMax(500); setMaxLipinskiViolations(1);
        setAdmetConfig(ADMET_PRESET_DEFAULTS.balanced); // ← NEW
    };


    useEffect(() => { pingBackend().then(setBackendAwake); }, []);

    useEffect(() => {
        if (pdbDebounceRef.current) clearTimeout(pdbDebounceRef.current);
        setPdbMeta(null);
        if (pdbId.length !== 4) return;
        pdbDebounceRef.current = setTimeout(async () => {
            setPdbMetaLoading(true);
            const meta = await fetchPdbMetadata(pdbId);
            setPdbMeta(meta);
            setPdbMetaLoading(false);
        }, 600);
    }, [pdbId]);

    const handleToxicityReportOnlyChange = (v: boolean) => setToxicityReportOnly(v);
    const handleDirectScoreOnlyChange = (v: boolean) => setDirectScoreOnly(v);

    const toxOnly = toxicityReportOnly && !directScoreOnly;
    const bothModes = toxicityReportOnly && directScoreOnly;
    const pdbRequired = !toxicityReportOnly;

    const handleSubmit = async (e: React.FormEvent) => {
        e.preventDefault();
        setError(null);
        if (!smiles.trim()) { setError("SMILES string is required."); return; }
        if (pdbRequired) {
            if (pdbMode === "id" && pdbId.trim().length !== 4) {
                setError("PDB ID must be exactly 4 characters (e.g. 1HSG)."); return;
            }
            if (pdbMode === "file" && !pdbFile) {
                setError("Please upload a .pdb file."); return;
            }
        }
        setSubmitting(true);
        try {
            const sharedOptions = {
                num_analogues: directScoreOnly ? 0 : numAnalogues,
                docking_speed: dockingSpeed,
                max_docking_compounds: maxDockingCompounds,
                binding_site_mode: bindingSiteMode,
                pipeline_steps: pipelineSteps,
                direct_score_only: directScoreOnly,
                mw_min: mwMin,
                mw_max: mwMax,
                max_lipinski_violations: maxLipinskiViolations,
                solubility_filter: solubilityFilter,
                toxicity_report_only: toxicityReportOnly,
                locked_scaffold_smarts: lockedScaffoldSmarts.trim() || undefined,
                admet_config: admetConfig, // ← NEW
            };
            let response;
            if (pdbMode === "file") {
                response = await submitJobWithFile(smiles, pdbFile, sharedOptions);
            } else {
                response = await submitJob({
                    smiles,
                    ...(pdbRequired && pdbId ? { pdb_id: pdbId.toUpperCase() } : {}),
                    ...sharedOptions,
                    binding_site_coords: bindingSiteMode === "coordinates" ? bindingSiteCoords : undefined,
                    binding_site_residues: bindingSiteMode === "residues" ? bindingSiteResidues : undefined,
                });
            }
            router.push(`/status/${response.job_id}`);
        } catch (err: unknown) {
            setError(err instanceof Error ? err.message : "Submission failed. Try again.");
            setSubmitting(false);
        }
    };

    const formatSavedTime = (d: Date): string => {
        const now = new Date();
        const diff = Math.round((now.getTime() - d.getTime()) / 1000);
        if (diff < 5) return "just now";
        if (diff < 60) return `${diff}s ago`;
        if (diff < 3600) return `${Math.floor(diff / 60)}m ago`;
        return d.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
    };

    return (
        <form onSubmit={handleSubmit} className="space-y-6">

            {memoryRestored && (
                <div className="flex items-center justify-between px-3 py-2 rounded-lg bg-gray-800/60 border border-gray-700/60 animate-slide-up">
                    <span className="flex items-center gap-2 text-xs text-gray-400">
                        <History className="w-3.5 h-3.5 text-emerald-500 flex-shrink-0" />
                        Previous settings restored — just update your SMILES and submit.
                    </span>
                    <button type="button" onClick={clearSavedSettings}
                        className="flex items-center gap-1 text-xs text-gray-600 hover:text-red-400 transition-colors ml-4 flex-shrink-0"
                        title="Reset all settings to defaults and clear saved memory">
                        <RotateCcw className="w-3 h-3" />Reset
                    </button>
                </div>
            )}

            {backendAwake === false && (
                <div className="flex items-start gap-3 p-4 rounded-xl bg-yellow-950/30 border border-yellow-800/50">
                    <AlertCircle className="w-4 h-4 text-yellow-400 mt-0.5 flex-shrink-0" />
                    <div className="text-sm">
                        <p className="text-yellow-300 font-medium">Backend is waking up</p>
                        <p className="text-yellow-600 mt-0.5">The Hugging Face Space may be sleeping. First request takes ~60s. You can still submit — it will queue automatically.</p>
                    </div>
                </div>
            )}

            {toxOnly && (
                <div className="flex items-start gap-3 p-4 rounded-xl bg-red-950/20 border border-red-800/40 animate-slide-up">
                    <TestTube className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" />
                    <div className="text-sm">
                        <p className="text-red-300 font-medium">Toxicity Report mode active</p>
                        <p className="text-red-700 mt-0.5">Only ADMET toxicity analysis will run. Docking and retrosynthesis are skipped. No protein structure required. Full report shown even if the compound fails ADMET.</p>
                    </div>
                </div>
            )}

            {bothModes && (
                <div className="flex items-start gap-3 p-4 rounded-xl bg-violet-950/20 border border-violet-800/40 animate-slide-up">
                    <Zap className="w-4 h-4 text-violet-400 mt-0.5 flex-shrink-0" />
                    <div className="text-sm">
                        <p className="text-violet-300 font-medium">Direct Score + Toxicity Report</p>
                        <p className="text-violet-700 mt-0.5">ADMET analysis on your base compound only — no analogues generated, no docking, no PDB required. Full toxicity profile is always written, even if the compound fails ADMET.</p>
                    </div>
                </div>
            )}

            <div className="card">
                <label className="block text-sm font-medium text-gray-300 mb-1">
                    Base Compound SMILES
                    <span className="ml-2 text-xs text-gray-500 font-normal">required</span>
                </label>
                <p className="text-xs text-gray-500 mb-3">
                    {bothModes
                        ? "This exact compound will be screened for ADMET toxicity (no analogues, no docking)."
                        : directScoreOnly
                            ? "This exact compound will be scored through the full pipeline (no analogues generated)."
                            : toxOnly
                                ? "This compound (and any analogues you request) will be screened for ADMET toxicity."
                                : "Enter the SMILES string of your base compound. Analogues will be generated from this scaffold."}
                </p>
                <div className="relative">
                    <textarea
                        className="input-base font-mono text-xs resize-none h-20"
                        value={smiles}
                        onChange={(e) => setSmiles(e.target.value)}
                        placeholder={MYRICETIN_SMILES}
                        spellCheck={false}
                    />
                    {smiles && (
                        <button type="button" onClick={() => setSmiles("")}
                            className="absolute top-2 right-2 text-gray-600 hover:text-gray-400 transition-colors">
                            <X className="w-3.5 h-3.5" />
                        </button>
                    )}
                </div>
                <div className="flex items-center gap-3 mt-2">
                    <button type="button" onClick={() => setSmiles(MYRICETIN_SMILES)}
                        className="text-xs text-emerald-500 hover:text-emerald-400 transition-colors">
                        Use myricetin (default)
                    </button>
                    <a href="https://pubchem.ncbi.nlm.nih.gov" target="_blank" rel="noopener noreferrer"
                        className="text-xs text-gray-500 hover:text-gray-400 flex items-center gap-1 transition-colors">
                        Find SMILES on PubChem <ExternalLink className="w-3 h-3" />
                    </a>
                </div>
            </div>

            {pdbRequired && (
                <div className="card">
                    <label className="block text-sm font-medium text-gray-300 mb-1">
                        Target Protein
                        <span className="ml-2 text-xs text-gray-500 font-normal">required</span>
                    </label>
                    <p className="text-xs text-gray-500 mb-3">Provide the protein structure to dock against.</p>

                    <div className="flex gap-1 p-1 bg-gray-800 rounded-lg w-fit mb-4">
                        {(["id", "file"] as const).map((mode) => (
                            <button key={mode} type="button" onClick={() => setPdbMode(mode)}
                                className={`px-4 py-1.5 rounded-md text-sm font-medium transition-all ${pdbMode === mode ? "bg-gray-700 text-gray-100" : "text-gray-500 hover:text-gray-300"}`}>
                                {mode === "id" ? "PDB ID" : "Upload File"}
                            </button>
                        ))}
                    </div>

                    {pdbMode === "id" && (
                        <div>
                            <div className="relative">
                                <input
                                    type="text"
                                    className="input-base font-mono uppercase pr-10"
                                    value={pdbId}
                                    onChange={(e) => setPdbId(e.target.value.toUpperCase().slice(0, 4))}
                                    placeholder="e.g. 1HSG"
                                    maxLength={4}
                                />
                                {pdbMetaLoading && (
                                    <Loader2 className="absolute right-3 top-2.5 w-4 h-4 text-gray-500 animate-spin" />
                                )}
                                {pdbMeta && !pdbMetaLoading && (
                                    <CheckCircle2 className="absolute right-3 top-2.5 w-4 h-4 text-emerald-500" />
                                )}
                            </div>
                            {pdbMeta && (
                                <div className="mt-2 p-3 rounded-lg bg-emerald-950/20 border border-emerald-900/40 animate-slide-up">
                                    <p className="text-sm text-emerald-300 font-medium truncate">{pdbMeta.title}</p>
                                    <div className="flex gap-4 mt-1">
                                        {pdbMeta.resolution_angstrom && (
                                            <span className="text-xs text-gray-500">
                                                Resolution <span className="text-gray-400">{pdbMeta.resolution_angstrom}Å</span>
                                            </span>
                                        )}
                                        <span className="text-xs text-gray-500">
                                            Chains <span className="text-gray-400">{pdbMeta.protein_chains}</span>
                                        </span>
                                        <a href={`https://www.rcsb.org/structure/${pdbId}`} target="_blank" rel="noopener noreferrer"
                                            className="text-xs text-emerald-600 hover:text-emerald-500 flex items-center gap-1">
                                            View on RCSB <ExternalLink className="w-3 h-3" />
                                        </a>
                                    </div>
                                </div>
                            )}
                            {pdbId.length === 4 && !pdbMeta && !pdbMetaLoading && (
                                <p className="mt-1 text-xs text-red-400">PDB ID not found. Check the ID at rcsb.org.</p>
                            )}
                        </div>
                    )}

                    {pdbMode === "file" && (
                        <div
                            onClick={() => fileInputRef.current?.click()}
                            className={`border-2 border-dashed rounded-xl p-6 text-center cursor-pointer transition-all ${pdbFile
                                ? "border-emerald-700 bg-emerald-950/20"
                                : "border-gray-700 hover:border-gray-600 bg-gray-800/20"
                                }`}>
                            <input
                                ref={fileInputRef}
                                type="file"
                                accept=".pdb"
                                className="hidden"
                                onChange={(e) => setPdbFile(e.target.files?.[0] ?? null)}
                            />
                            {pdbFile ? (
                                <div className="flex items-center justify-center gap-2">
                                    <CheckCircle2 className="w-5 h-5 text-emerald-500" />
                                    <div className="text-left">
                                        <p className="text-sm text-emerald-300 font-medium">{pdbFile.name}</p>
                                        <p className="text-xs text-gray-500">{(pdbFile.size / 1024).toFixed(1)} KB</p>
                                    </div>
                                    <button type="button"
                                        onClick={(e) => { e.stopPropagation(); setPdbFile(null); }}
                                        className="ml-2 text-gray-600 hover:text-gray-400">
                                        <X className="w-4 h-4" />
                                    </button>
                                </div>
                            ) : (
                                <>
                                    <Upload className="w-8 h-8 text-gray-600 mx-auto mb-2" />
                                    <p className="text-sm text-gray-400">Click to upload a .pdb file</p>
                                    <p className="text-xs text-gray-600 mt-1">Download from rcsb.org → Format: PDB → Download</p>
                                </>
                            )}
                        </div>
                    )}
                </div>
            )}

            {/* AdvancedSettings now receives admetConfig + setter */}
            <AdvancedSettings
                numAnalogues={numAnalogues}
                onNumAnaloguesChange={setNumAnalogues}
                directScoreOnly={directScoreOnly}
                onDirectScoreOnlyChange={handleDirectScoreOnlyChange}
                toxicityReportOnly={toxicityReportOnly}
                onToxicityReportOnlyChange={handleToxicityReportOnlyChange}
                solubilityFilter={solubilityFilter}
                onSolubilityFilterChange={setSolubilityFilter}
                lockedScaffoldSmarts={lockedScaffoldSmarts}
                onLockedScaffoldSmartsChange={setLockedScaffoldSmarts}
                pipelineSteps={pipelineSteps}
                onPipelineStepsChange={setPipelineSteps}
                dockingSpeed={dockingSpeed}
                onDockingSpeedChange={setDockingSpeed}
                maxDockingCompounds={maxDockingCompounds}
                onMaxDockingCompoundsChange={setMaxDockingCompounds}
                bindingSiteMode={bindingSiteMode}
                onBindingSiteModeChange={setBindingSiteMode}
                bindingSiteCoords={bindingSiteCoords}
                onBindingSiteCoordsChange={setBindingSiteCoords}
                bindingSiteResidues={bindingSiteResidues}
                onBindingSiteResiduesChange={setBindingSiteResidues}
                mwMin={mwMin}
                mwMax={mwMax}
                onMwMinChange={setMwMin}
                onMwMaxChange={setMwMax}
                maxLipinskiViolations={maxLipinskiViolations}
                onMaxLipinskiViolationsChange={setMaxLipinskiViolations}
                admetConfig={admetConfig}
                onAdmetConfigChange={setAdmetConfig}
            />

            {error && (
                <div className="flex items-start gap-3 p-4 rounded-xl bg-red-950/30 border border-red-800/50 animate-slide-up">
                    <AlertCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" />
                    <p className="text-sm text-red-300">{error}</p>
                </div>
            )}

            <button type="submit" disabled={submitting} className="btn-primary w-full py-3 text-base">
                {submitting
                    ? <><Loader2 className="w-5 h-5 animate-spin" /> Submitting pipeline job...</>
                    : bothModes
                        ? <><Zap className="w-5 h-5" /> Run Toxicity Report (Single Compound)</>
                        : toxOnly
                            ? <><TestTube className="w-5 h-5" /> Run Toxicity Report</>
                            : directScoreOnly
                                ? <><FlaskConical className="w-5 h-5" /> Score This Compound Directly</>
                                : <><FlaskConical className="w-5 h-5" /> Run CADD Pipeline</>
                }
            </button>

            {!submitting && (
                <div className="flex items-center justify-between">
                    <p className="text-xs text-gray-600">
                        {bothModes
                            ? "ADMET toxicity report on base compound only — no PDB needed, ~2–5 minutes."
                            : toxOnly
                                ? "Toxicity report only — ADMET screening completes in ~2–5 minutes."
                                : directScoreOnly
                                    ? "Direct scoring skips analogue generation — results in 5–15 minutes."
                                    : numAnalogues >= 500
                                        ? `${numAnalogues} analogues requested — pipeline may take 60–90 minutes.`
                                        : "Pipeline takes 30–60 minutes depending on settings. You'll get a tracking URL immediately."}
                    </p>
                    {savedAt && (
                        <span className="flex items-center gap-1.5 text-xs text-gray-600 ml-4 flex-shrink-0 whitespace-nowrap">
                            <History className="w-3 h-3 text-emerald-600" />
                            Saved {formatSavedTime(savedAt)}
                        </span>
                    )}
                </div>
            )}

        </form>
    );
}
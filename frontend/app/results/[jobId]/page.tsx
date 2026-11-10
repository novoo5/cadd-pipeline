"use client";

import { useEffect, useState } from "react";
import { useParams, useRouter } from "next/navigation";
import { Home, ArrowLeft, Loader2, MapPin, AlertCircle } from "lucide-react";
import ResultsTable from "@/components/ResultsTable";
import { getJobResults, type JobResultsResponse } from "@/lib/api";

export default function ResultsPage() {
    const { jobId } = useParams<{ jobId: string }>();
    const router = useRouter();

    const [results, setResults] = useState<JobResultsResponse | null>(null);
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        const load = async () => {
            try {
                const data = await getJobResults(jobId);
                if (data === null) {
                    router.replace(`/status/${jobId}`);
                    return;
                }
                setResults(data);
            } catch (e: unknown) {
                setError(e instanceof Error ? e.message : "Failed to load results.");
            } finally {
                setLoading(false);
            }
        };
        load();
    }, [jobId, router]);

    if (loading) {
        return (
            <div className="flex flex-col items-center justify-center min-h-64 gap-4">
                <Loader2 className="w-8 h-8 text-emerald-500 animate-spin" />
                <p className="text-sm text-gray-500">Loading results...</p>
            </div>
        );
    }

    if (error) {
        return (
            <div className="max-w-2xl mx-auto">
                <div className="card">
                    <p className="text-red-400 text-sm mb-4">{error}</p>
                    <div className="flex gap-2">
                        <button onClick={() => router.push(`/status/${jobId}`)} className="btn-secondary text-sm">
                            <ArrowLeft className="w-4 h-4" /> Back to Status
                        </button>
                        <button onClick={() => router.push("/")} className="btn-secondary text-sm">
                            <Home className="w-4 h-4" /> New Job
                        </button>
                    </div>
                </div>
            </div>
        );
    }

    if (!results) return null;

    const bsi = (results as any).binding_site_info;

    const detectionBadge = bsi ? {
        native_ligand: { label: `Native ligand: ${bsi.detected_ligand_name}`, color: "text-emerald-400 bg-emerald-950/40 border-emerald-800/40" },
        protein_centroid: { label: "Fallback: protein centroid (no native ligand found)", color: "text-yellow-400 bg-yellow-950/40 border-yellow-800/40" },
        user_coordinates: { label: "User-supplied coordinates", color: "text-blue-400 bg-blue-950/40 border-blue-800/40" },
        user_residues: { label: `Residue range: ${bsi.detected_ligand_name}`, color: "text-blue-400 bg-blue-950/40 border-blue-800/40" },
    }[bsi.detection_mode as string] : null;

    return (
        <div className="space-y-6">

            {/* ── Header ──────────────────────────────────────────────────────── */}
            <div>
                <div className="flex items-center gap-2 text-xs text-gray-600 mb-3">
                    <a href="/" className="hover:text-gray-400 transition-colors flex items-center gap-1">
                        <Home className="w-3 h-3" /> Home
                    </a>
                    <span>/</span>
                    <a href={`/status/${jobId}`} className="hover:text-gray-400 transition-colors">Status</a>
                    <span>/</span>
                    <span className="text-gray-500">Results</span>
                </div>

                <div className="flex items-start justify-between gap-4">
                    <div>
                        <h1 className="text-2xl font-bold text-gray-100">Pipeline Results</h1>
                        <p className="text-sm text-gray-500 mt-1">
                            {results.final_ranked_compounds.length} compounds ranked ·{" "}
                            <span className="font-mono text-xs text-gray-600">
                                Base: {results.base_smiles.slice(0, 30)}
                                {results.base_smiles.length > 30 ? "..." : ""}
                            </span>
                        </p>
                    </div>
                    <div className="flex gap-2 flex-shrink-0">
                        <button onClick={() => router.push(`/status/${jobId}`)} className="btn-secondary text-xs py-1.5">
                            <ArrowLeft className="w-3.5 h-3.5" /> Status
                        </button>
                        <button onClick={() => router.push("/")} className="btn-secondary text-xs py-1.5">
                            <Home className="w-3.5 h-3.5" /> New Job
                        </button>
                    </div>
                </div>

                <div className="flex gap-4 mt-3 text-xs text-gray-600">
                    <span>Started: {new Date(results.created_at).toLocaleString()}</span>
                    <span>·</span>
                    <span>Completed: {new Date(results.completed_at).toLocaleString()}</span>
                    <span>·</span>
                    <span>
                        Duration:{" "}
                        {Math.round(
                            (new Date(results.completed_at).getTime() -
                                new Date(results.created_at).getTime()) / 60000
                        )} min
                    </span>
                </div>
            </div>

            {/* ── Binding Site Info Card ───────────────────────────────────────── */}
            {bsi && (
                <div className={`p-4 rounded-xl border ${detectionBadge?.color ?? "border-gray-700 bg-gray-800/40"}`}>
                    <div className="flex items-center gap-2 mb-3">
                        <MapPin className="w-4 h-4 flex-shrink-0" />
                        <p className="text-sm font-medium text-gray-200">Docking Grid Center</p>
                        <span className={`text-xs px-2 py-0.5 rounded-full border ${detectionBadge?.color}`}>
                            {detectionBadge?.label}
                        </span>
                    </div>

                    {bsi.detection_mode === "protein_centroid" && (
                        <div className="flex items-start gap-2 mb-3 p-2.5 rounded-lg bg-yellow-950/30 border border-yellow-800/30">
                            <AlertCircle className="w-3.5 h-3.5 text-yellow-400 flex-shrink-0 mt-0.5" />
                            <p className="text-xs text-yellow-300/80 leading-relaxed">
                                No co-crystallized ligand was found in this PDB file. The docking grid was placed at the
                                geometric center of the full protein — this is <strong>blind docking</strong> and may be less
                                accurate. For best results, use a PDB with a native ligand (e.g. 4COX has Ibuprofen, 2HU4 has Tamiflu).
                            </p>
                        </div>
                    )}

                    {bsi.detection_mode === "native_ligand" && (
                        <p className="text-xs text-gray-500 mb-3 leading-relaxed">
                            The pipeline automatically detected the co-crystallized native ligand <strong className="text-gray-300">{bsi.detected_ligand_name}</strong> in the
                            PDB file and centered the docking grid on its binding pocket. All analogues were docked into
                            this exact pocket. You can verify or override these coordinates below.
                        </p>
                    )}

                    <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
                        {[
                            { label: "Center X", value: bsi.center_x.toFixed(3), unit: "Å" },
                            { label: "Center Y", value: bsi.center_y.toFixed(3), unit: "Å" },
                            { label: "Center Z", value: bsi.center_z.toFixed(3), unit: "Å" },
                            { label: "Box Size", value: bsi.box_size.toFixed(1), unit: "Å" },
                        ].map(({ label, value, unit }) => (
                            <div key={label} className="p-2.5 rounded-lg bg-gray-900/60 border border-gray-700/50 text-center">
                                <p className="text-xs text-gray-500 mb-0.5">{label}</p>
                                <p className="font-mono text-sm font-medium text-gray-200">
                                    {value} <span className="text-gray-500 text-xs">{unit}</span>
                                </p>
                            </div>
                        ))}
                    </div>

                    <p className="text-xs text-gray-600 mt-2.5">
                        To use different coordinates next time, select <strong className="text-gray-500">Manual Coordinates</strong> mode
                        in the job form and enter the values above as a starting point.
                    </p>
                </div>
            )}

            {/* ── Results table ────────────────────────────────────────────────── */}
            <ResultsTable results={results} />
        </div>
    );
}

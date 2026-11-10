"use client";

import { useEffect, useState, useCallback } from "react";
import { useParams, useRouter } from "next/navigation";
import { ArrowRight, RefreshCw, Home } from "lucide-react";
import StatusTracker from "@/components/StatusTracker";
import { getJobStatus, type JobStatusResponse } from "@/lib/api";

const POLL_INTERVAL_MS = 8000;

export default function StatusPage() {
    const { jobId } = useParams<{ jobId: string }>();
    const router = useRouter();

    const [status, setStatus] = useState<JobStatusResponse | null>(null);
    const [error, setError] = useState<string | null>(null);
    const [lastPoll, setLastPoll] = useState<Date | null>(null);

    const poll = useCallback(async () => {
        try {
            const data = await getJobStatus(jobId);
            setStatus(data);
            setLastPoll(new Date());
            setError(null);

            // Auto-navigate to results when done
            if (data.status === "done") {
                setTimeout(() => router.push(`/results/${jobId}`), 1500);
            }
        } catch (e: unknown) {
            setError(e instanceof Error ? e.message : "Failed to fetch status.");
        }
    }, [jobId, router]);

    // Initial fetch + polling interval
    useEffect(() => {
        poll();
        const interval = setInterval(poll, POLL_INTERVAL_MS);
        return () => clearInterval(interval);
    }, [poll]);

    return (
        <div className="max-w-2xl mx-auto space-y-6">

            {/* ── Header ─────────────────────────────────────────────────────── */}
            <div>
                <div className="flex items-center gap-2 text-xs text-gray-600 mb-3">
                    <a href="/" className="hover:text-gray-400 transition-colors flex items-center gap-1">
                        <Home className="w-3 h-3" /> Home
                    </a>
                    <span>/</span>
                    <span className="text-gray-500">Job Status</span>
                </div>
                <h1 className="text-2xl font-bold text-gray-100">Pipeline Status</h1>
                <p className="text-sm text-gray-500 mt-1 font-mono">
                    Job ID: <span className="text-gray-400">{jobId}</span>
                </p>
            </div>

            {/* ── Status tracker ────────────────────────────────────────────── */}
            {status ? (
                <StatusTracker
                    steps={status.steps}
                    overallStatus={status.status}
                    error={status.error}
                />
            ) : error ? (
                <div className="card">
                    <p className="text-sm text-red-400">{error}</p>
                </div>
            ) : (
                // Loading skeleton
                <div className="space-y-2">
                    {Array.from({ length: 7 }).map((_, i) => (
                        <div key={i} className="skeleton h-14 rounded-xl" />
                    ))}
                </div>
            )}

            {/* ── Actions ───────────────────────────────────────────────────── */}
            <div className="flex items-center justify-between pt-2">
                <div className="flex items-center gap-2 text-xs text-gray-600">
                    {lastPoll && (
                        <>
                            <RefreshCw className="w-3 h-3" />
                            Last updated {lastPoll.toLocaleTimeString()} · auto-refreshes every 8s
                        </>
                    )}
                </div>

                <div className="flex gap-2">
                    <button onClick={() => router.push("/")} className="btn-secondary text-xs py-1.5">
                        <Home className="w-3.5 h-3.5" />
                        New Job
                    </button>
                    {status?.status === "done" && (
                        <button
                            onClick={() => router.push(`/results/${jobId}`)}
                            className="btn-primary text-xs py-1.5"
                        >
                            View Results
                            <ArrowRight className="w-3.5 h-3.5" />
                        </button>
                    )}
                </div>
            </div>

            {/* ── Job info card ─────────────────────────────────────────────── */}
            {status && (
                <div className="card text-xs text-gray-600 space-y-1">
                    <div className="flex justify-between">
                        <span>Submitted</span>
                        <span className="text-gray-500">
                            {new Date(status.created_at).toLocaleString()}
                        </span>
                    </div>
                    <div className="flex justify-between">
                        <span>Last updated</span>
                        <span className="text-gray-500">
                            {new Date(status.updated_at).toLocaleString()}
                        </span>
                    </div>
                    <div className="flex justify-between">
                        <span>Bookmark this page to check back later</span>
                        <span className="text-emerald-600 font-medium">↑ save the URL</span>
                    </div>
                </div>
            )}
        </div>
    );
}

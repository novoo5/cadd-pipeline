"use client";

import { CheckCircle2, Circle, Loader2, XCircle, MinusCircle, Clock } from "lucide-react";
import type { PipelineStepInfo, StepStatus } from "@/lib/api";
import { formatDuration } from "@/lib/api";

interface StatusTrackerProps {
    steps: PipelineStepInfo[];
    overallStatus: "queued" | "running" | "done" | "failed";
    error?: string;
}

const STEP_ICONS: Record<StepStatus, React.ReactNode> = {
    waiting: <Circle className="w-4 h-4 text-gray-600" />,
    running: <Loader2 className="w-4 h-4 text-blue-400 animate-spin" />,
    done: <CheckCircle2 className="w-4 h-4 text-emerald-400" />,
    skipped: <MinusCircle className="w-4 h-4 text-gray-600" />,
    failed: <XCircle className="w-4 h-4 text-red-400" />,
};

const STEP_CLASS: Record<StepStatus, string> = {
    waiting: "step-waiting",
    running: "step-running",
    done: "step-done",
    skipped: "step-skipped",
    failed: "step-failed",
};

const OVERALL_LABEL: Record<string, { text: string; color: string }> = {
    queued: { text: "Queued — waiting for worker", color: "text-gray-400" },
    running: { text: "Running pipeline...", color: "text-blue-400" },
    done: { text: "Pipeline complete ✓", color: "text-emerald-400" },
    failed: { text: "Pipeline failed", color: "text-red-400" },
};

export default function StatusTracker({
    steps,
    overallStatus,
    error,
}: StatusTrackerProps) {
    const doneCount = steps.filter((s) => s.status === "done").length;
    const totalActive = steps.filter((s) => s.status !== "skipped").length;
    const progressPct = totalActive > 0 ? (doneCount / totalActive) * 100 : 0;
    const label = OVERALL_LABEL[overallStatus] ?? OVERALL_LABEL.queued;

    return (
        <div className="space-y-4">

            {/* ── Overall status header ──────────────────────────────────────── */}
            <div className="card">
                <div className="flex items-center justify-between mb-3">
                    <div className="flex items-center gap-2">
                        {overallStatus === "running" && (
                            <Loader2 className="w-4 h-4 text-blue-400 animate-spin" />
                        )}
                        {overallStatus === "done" && (
                            <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                        )}
                        {overallStatus === "failed" && (
                            <XCircle className="w-4 h-4 text-red-400" />
                        )}
                        {overallStatus === "queued" && (
                            <Clock className="w-4 h-4 text-gray-500" />
                        )}
                        <span className={`text-sm font-medium ${label.color}`}>
                            {label.text}
                        </span>
                    </div>
                    <span className="text-xs text-gray-600">
                        {doneCount}/{totalActive} steps done
                    </span>
                </div>

                {/* Overall pipeline progress bar */}
                <div className="h-1.5 bg-gray-800 rounded-full overflow-hidden">
                    <div
                        className={`h-full rounded-full transition-all duration-700 ease-out ${overallStatus === "failed"
                                ? "bg-red-600"
                                : overallStatus === "done"
                                    ? "bg-emerald-500"
                                    : "bg-blue-500"
                            }`}
                        style={{ width: `${overallStatus === "done" ? 100 : progressPct}%` }}
                    />
                </div>
            </div>

            {/* ── Step list ─────────────────────────────────────────────────── */}
            <div className="space-y-2">
                {steps.map((step, i) => {
                    // ← derive step-level progress percentage if available
                    const hasStepProgress =
                        step.status === "running" &&
                        step.progress_total != null &&
                        step.progress_total > 0;
                    const stepPct = hasStepProgress
                        ? Math.round((step.progress_current! / step.progress_total!) * 100)
                        : 0;

                    return (
                        <div
                            key={i}
                            className={`flex items-start gap-3 p-3.5 rounded-xl border transition-all animate-fade-in ${STEP_CLASS[step.status]}`}
                        >
                            {/* Step number + icon */}
                            <div className="flex items-center gap-2 flex-shrink-0 pt-0.5">
                                <span className="text-xs text-gray-700 font-mono w-4 text-right">
                                    {i + 1}
                                </span>
                                {STEP_ICONS[step.status]}
                            </div>

                            {/* Step content */}
                            <div className="flex-1 min-w-0">
                                <div className="flex items-center justify-between gap-2">
                                    <span className="text-sm font-medium">{step.name}</span>
                                    {step.duration_seconds != null ? (
                                        <span className="text-xs text-gray-600 flex-shrink-0 flex items-center gap-1">
                                            <Clock className="w-3 h-3" />
                                            {formatDuration(step.duration_seconds)}
                                        </span>
                                    ) : hasStepProgress ? (
                                        // ← live counter in top-right while running
                                        <span className="text-xs text-blue-400 flex-shrink-0 font-mono tabular-nums">
                                            {step.progress_current}/{step.progress_total}
                                        </span>
                                    ) : null}
                                </div>

                                {step.message && (
                                    <p className="text-xs mt-0.5 opacity-80 truncate">
                                        {step.message}
                                    </p>
                                )}
                                {step.status === "running" && !step.message && !hasStepProgress && (
                                    <p className="text-xs mt-0.5 opacity-60">Processing...</p>
                                )}

                                {/* ← Step-level progress bar (only while running + progress exists) */}
                                {hasStepProgress && (
                                    <div className="mt-2">
                                        <div className="flex justify-between items-center text-[10px] text-gray-500 mb-1">
                                            <span className="opacity-70">
                                                {step.name.includes("Docking")
                                                    ? "Compounds docked"
                                                    : "Analogues generated"}
                                            </span>
                                            <span className="font-mono tabular-nums">{stepPct}%</span>
                                        </div>
                                        <div className="w-full h-1.5 bg-gray-800 rounded-full overflow-hidden">
                                            <div
                                                className="h-full bg-blue-500 rounded-full transition-all duration-500 ease-out"
                                                style={{ width: `${stepPct}%` }}
                                            />
                                        </div>
                                    </div>
                                )}
                            </div>
                        </div>
                    );
                })}
            </div>

            {/* ── Error detail ──────────────────────────────────────────────── */}
            {error && (
                <div className="p-4 rounded-xl bg-red-950/30 border border-red-800/50 animate-slide-up">
                    <p className="text-xs font-medium text-red-400 mb-1">Error Detail</p>
                    <p className="text-xs text-red-300 font-mono break-all">{error}</p>
                </div>
            )}
        </div>
    );
}

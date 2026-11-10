'use client'

import { useState, useMemo, useCallback } from 'react'
import {
    ChevronDown, ChevronUp, Copy, CheckCheck, FlaskConical, Activity,
    Dna, GitBranch, Info, AlertTriangle, AlertCircle, Lightbulb,
    Download, Loader2, BookOpen, Layers, Shield, CheckCircle2, Clock,
    Thermometer, ArrowRight, FileText, Percent, TestTube,
    SlidersHorizontal, XCircle, Filter,
} from 'lucide-react'
import type { CompoundResult, ScoreBreakdown, ADMETFlagDetail, DockingFileType } from '@/lib/api'
import {
    getAffinityColor, formatProbability,
    getScoreBreakdown, getFlagSeverityColor, downloadDockingFile,
} from '@/lib/api'

// ── Extended ADMET types ──────────────────────────────────────────────────────

interface ADMETNarrativeBlock {
    title: string
    level: 'info' | 'good' | 'warning' | 'danger'
    body: string
}

interface ADMETEndpointValue {
    key: string
    label: string
    category: 'absorption' | 'distribution' | 'metabolism' | 'excretion' | 'toxicity' | 'physchem' | 'other'
    value: number
    display_value: string
    unit?: string | null
    interpretation?: string | null
    threshold_applied?: string | null
    severity?: string | null
    hard_fail: boolean
    triggered: boolean
}

interface ADMETThresholdConfig {
    enabled: boolean
    cutoff: number
    direction: 'above' | 'below'
    severity_high?: number | null
    hard_fail: boolean
    threshold_str?: string | null
    implication?: string | null
    recommendation?: string | null
}

type OverallRisk = 'low' | 'moderate' | 'high'

interface ExtendedADMET {
    passed: boolean
    overall_risk?: OverallRisk
    herg_inhibition: number
    caco2_permeability: number
    bbb_penetration: number
    hepatotoxicity: number
    oral_bioavailability: number
    flags: ADMETFlagDetail[]
    flag_summary: string[]
    decision_basis?: string[]
    narrative?: ADMETNarrativeBlock[]
    endpoint_table?: ADMETEndpointValue[]
    extra_properties?: Record<string, number>
    thresholds_used?: Record<string, ADMETThresholdConfig>
}

interface ExtendedRetrosynthesisStep {
    step_number: number
    reaction_smarts: string
    reaction_name?: string | null
    reaction_type?: string | null
    starting_materials: string[]
    product_smiles?: string | null
    schematic?: string | null
    reagents?: string[]
    solvents?: string[]
    conditions?: string | null
    temperature?: string | null
    duration?: string | null
    yield_estimate?: string | null
    protocol_text?: string | null
    confidence: number
}

interface ExtendedRetrosynthesisResult {
    feasible: boolean
    num_steps: number
    sa_score?: number
    route: ExtendedRetrosynthesisStep[]
    complexity_score: number
    difficulty_label?: string
    synthesis_summary?: string | null
    estimated_total_yield?: string | null
}

// ── Breakdown item shape from backend v4 ─────────────────────────────────────
// Each scored entry has: score_0_to_1, raw, contribution, max_possible,
// base_weight, effective_weight, step_active
// binding_prefilter has gate_only:true, contribution:0, max_possible:0
// weight_redistribution and penalties are info-only objects
// final_score is a plain number

interface BreakdownEntry {
    score_0_to_1?: number
    raw?: string
    contribution?: number
    max_possible?: number        // = effective_weight (whole-number pts)
    effective_weight?: number
    base_weight?: number
    step_active?: boolean
    explanation?: string
    gate_only?: boolean          // binding_prefilter sentinel
    passed?: boolean             // binding_prefilter gate result
}

// ── Constants ─────────────────────────────────────────────────────────────────

const CATEGORY_ORDER = ['absorption', 'distribution', 'metabolism', 'excretion', 'toxicity', 'physchem', 'other'] as const

const CATEGORY_LABELS: Record<string, string> = {
    absorption: 'Absorption', distribution: 'Distribution', metabolism: 'Metabolism',
    excretion: 'Excretion', toxicity: 'Toxicity', physchem: 'Phys. Chem.', other: 'Other',
}

const THRESHOLD_KEY_LABELS: Record<string, string> = {
    herg_inhibition: 'hERG Inhibition', hepatotoxicity: 'Hepatotoxicity',
    caco2_permeability: 'Caco-2 Permeability', oral_bioavailability: 'Oral Bioavailability',
    bbb_penetration: 'BBB Penetration',
}

const BREAKDOWN_TIPS: Record<string, string> = {
    docking_affinity: 'Up to 51 pts. Normalised from −3 (no binding) to −12 kcal/mol (exceptional). More negative affinity = more points.',
    admet_safety: 'Up to 28 pts. 0 flags = full score. Each toxicity flag reduces the score. hERG and hepatotoxicity apply a ×2 hard-fail multiplier.',
    drug_likeness: 'Up to 11 pts. Weighted combination: Lipinski violations (50%) + logS solubility (30%) + logP optimality (20%). LogS ≥ −3: Excellent; ≥ −4: Good; ≥ −5: Borderline; < −6: Insoluble.',
    synthesis_ease: 'Up to 10 pts. 60% normalised SA-Score complexity + 40% step-count penalty.',
}

// Keys rendered as informational panels, not scored rows
const INFO_ONLY_KEYS = new Set(['weight_redistribution', 'penalties', 'final_score'])
// binding_prefilter is handled separately (gate_only block)

const DIFFICULTY_STYLES: Record<string, string> = {
    easy: 'text-emerald-400', moderate: 'text-yellow-400',
    hard: 'text-orange-400', infeasible: 'text-red-400', unknown: 'text-gray-500',
}

const REACTION_TYPE_LABELS: Record<string, string> = {
    transition_metal_catalysis: 'Transition Metal Catalysis', condensation: 'Condensation',
    condensation_reduction: 'Condensation + Reduction', nucleophilic_substitution: 'Nucleophilic Substitution',
    cyclisation: 'Cyclisation', multicomponent: 'Multicomponent',
    heteroatom_coupling: 'Heteroatom Coupling', carbon_coupling: 'C–C Coupling',
    fgi: 'Functional Group Interconversion',
}

// Display order for scored rows
const SCORED_KEY_ORDER = ['docking_affinity', 'admet_safety', 'drug_likeness', 'synthesis_ease']

// ── Style helpers ─────────────────────────────────────────────────────────────

const getRiskStyle = (risk: OverallRisk) => {
    switch (risk) {
        case 'high': return 'bg-red-950/50 border-red-800 text-red-400'
        case 'moderate': return 'bg-yellow-950/50 border-yellow-800 text-yellow-400'
        default: return 'bg-emerald-950/50 border-emerald-800 text-emerald-400'
    }
}

const getNarrativeStyle = (level: ADMETNarrativeBlock['level']) => {
    switch (level) {
        case 'danger': return { wrap: 'bg-red-950/30 border-red-900/60', title: 'text-red-300', body: 'text-red-200/70' }
        case 'warning': return { wrap: 'bg-yellow-950/30 border-yellow-900/60', title: 'text-yellow-300', body: 'text-yellow-200/70' }
        case 'good': return { wrap: 'bg-emerald-950/30 border-emerald-900/60', title: 'text-emerald-300', body: 'text-emerald-200/70' }
        default: return { wrap: 'bg-gray-800/50 border-gray-700', title: 'text-gray-300', body: 'text-gray-400' }
    }
}

const getNarrativeIcon = (level: ADMETNarrativeBlock['level']) => {
    switch (level) {
        case 'danger': return <AlertCircle className="w-3.5 h-3.5 flex-shrink-0 text-red-500 mt-0.5" />
        case 'warning': return <AlertTriangle className="w-3.5 h-3.5 flex-shrink-0 text-yellow-500 mt-0.5" />
        case 'good': return <CheckCircle2 className="w-3.5 h-3.5 flex-shrink-0 text-emerald-500 mt-0.5" />
        default: return <Info className="w-3.5 h-3.5 flex-shrink-0 text-gray-500 mt-0.5" />
    }
}

const getEndpointStatusClass = (ep: ADMETEndpointValue) => {
    if (!ep.triggered) return 'text-emerald-400'
    if (ep.severity === 'high') return 'text-red-400'
    return 'text-yellow-400'
}

const getAffinityLabel = (kcal: number): { label: string; color: string } => {
    if (kcal <= -9.0) return { label: 'Outstanding', color: 'text-emerald-300' }
    if (kcal <= -7.0) return { label: 'Strong', color: 'text-emerald-400' }
    if (kcal <= -5.0) return { label: 'Moderate', color: 'text-yellow-400' }
    return { label: 'Weak', color: 'text-red-400' }
}

const getComplexityLabel = (score: number): { label: string; color: string } => {
    if (score <= 15) return { label: 'Very Easy', color: 'text-emerald-400' }
    if (score <= 25) return { label: 'Moderate', color: 'text-yellow-400' }
    if (score <= 40) return { label: 'Difficult', color: 'text-orange-400' }
    return { label: 'Infeasible', color: 'text-red-400' }
}

// ── Tooltip ───────────────────────────────────────────────────────────────────

interface TooltipBubbleProps {
    text: string
    widthClass?: string
    align?: 'center' | 'left'
    iconClassName?: string
    className?: string
}

function TooltipBubble({
    text,
    widthClass = 'w-56',
    align = 'center',
    iconClassName = 'w-3 h-3 text-gray-600 hover:text-gray-400',
    className = '',
}: TooltipBubbleProps) {
    const [open, setOpen] = useState(false)
    const bubblePosition = align === 'left' ? 'left-0' : 'left-1/2 -translate-x-1/2'

    return (
        <span
            className={`relative inline-flex items-center ml-1 ${className}`}
            onMouseEnter={() => setOpen(true)}
            onMouseLeave={() => setOpen(false)}
            onFocus={() => setOpen(true)}
            onBlur={() => setOpen(false)}
        >
            <button
                type="button"
                className="inline-flex items-center justify-center cursor-help"
                onClick={(e) => { e.stopPropagation(); setOpen(v => !v) }}
                aria-label="Show info"
            >
                <Info className={`${iconClassName} transition-colors`} />
            </button>
            {open && (
                <span className={`absolute bottom-full ${bubblePosition} mb-1.5 ${widthClass} rounded-lg border border-gray-700 bg-gray-900 px-2.5 py-2 text-xs leading-relaxed text-gray-300 shadow-xl z-50 pointer-events-none`}>
                    {text}
                </span>
            )}
        </span>
    )
}

function Tip({ text }: { text: string }) {
    return <TooltipBubble text={text} widthClass="w-56" align="center" />
}

// ── Drug Likeness Chips ───────────────────────────────────────────────────────

interface DrugLikenessChip {
    label: string
    value: string
    status: 'good' | 'warn' | 'bad' | 'neutral'
}

function parseDrugLikenessRaw(raw: string): DrugLikenessChip[] {
    const chips: DrugLikenessChip[] = []
    const parts = raw.split('|').map(s => s.trim()).filter(Boolean)

    for (const part of parts) {
        const vMatch = part.match(/^Violations=(\d+)$/)
        if (vMatch) {
            const n = parseInt(vMatch[1])
            chips.push({ label: 'Violations', value: vMatch[1], status: n === 0 ? 'good' : n === 1 ? 'warn' : 'bad' })
            continue
        }

        const logsMatch = part.match(/^LogS=([-\d.]+)(?:\s*\(([^)]+)\))?$/)
        if (logsMatch) {
            const val = parseFloat(logsMatch[1])
            const label = logsMatch[2] ?? ''
            const lowerLabel = label.toLowerCase()
            let status: DrugLikenessChip['status'] = 'neutral'
            // v4 thresholds: ≥-3 Excellent, ≥-4 Good, ≥-5 Borderline, ≥-6 Poor, <-6 Insoluble
            if (lowerLabel === 'high' || val >= -3) status = 'good'
            else if (val >= -4) status = 'good'
            else if (val >= -5) status = 'warn'
            else status = 'bad'
            chips.push({ label: 'LogS', value: `${logsMatch[1]}${label ? ` (${label})` : ''}`, status })
            continue
        }

        const logpMatch = part.match(/^LogP=([-\d.]+)$/)
        if (logpMatch) {
            const val = parseFloat(logpMatch[1])
            // ideal 1.0-3.0, ok 0-5, penalty above 5
            chips.push({
                label: 'LogP',
                value: logpMatch[1],
                status: val >= 1.0 && val <= 3.0 ? 'good' : val >= 0 && val <= 5.0 ? 'warn' : 'bad',
            })
            continue
        }

        const mwMatch = part.match(/^MW=([\d.]+)\s*(Da)?$/)
        if (mwMatch) {
            const val = parseFloat(mwMatch[1])
            chips.push({ label: 'MW', value: `${mwMatch[1]} Da`, status: val <= 500 ? 'good' : val <= 600 ? 'warn' : 'bad' })
            continue
        }

        chips.push({ label: part, value: '', status: 'neutral' })
    }

    return chips
}

const CHIP_STYLES: Record<DrugLikenessChip['status'], string> = {
    good: 'bg-emerald-950/50 border-emerald-800/60 text-emerald-400',
    warn: 'bg-yellow-950/40 border-yellow-800/50 text-yellow-400',
    bad: 'bg-red-950/40 border-red-800/50 text-red-400',
    neutral: 'bg-gray-800/50 border-gray-700 text-gray-400',
}

function DrugLikenessChips({ raw }: { raw: string }) {
    if (!raw) return null
    const chips = parseDrugLikenessRaw(raw)
    if (chips.length === 0) return <p className="text-[10px] text-gray-600 leading-relaxed">{raw}</p>
    return (
        <div className="flex flex-wrap gap-1.5 mt-1">
            {chips.map((chip, i) => (
                <span key={i} className={`inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-[10px] font-medium border ${CHIP_STYLES[chip.status]}`}>
                    {chip.status === 'good' && <CheckCircle2 className="w-2.5 h-2.5 flex-shrink-0" />}
                    {chip.status === 'warn' && <AlertTriangle className="w-2.5 h-2.5 flex-shrink-0" />}
                    {chip.status === 'bad' && <AlertCircle className="w-2.5 h-2.5 flex-shrink-0" />}
                    <span className="text-gray-500 mr-0.5">{chip.label}</span>
                    {chip.value && <span className="font-mono">{chip.value}</span>}
                </span>
            ))}
        </div>
    )
}

// ── ADMET Flag Card ───────────────────────────────────────────────────────────

function ADMETFlagCard({ flag }: { flag: ADMETFlagDetail }) {
    const [open, setOpen] = useState(false)
    const severityClasses = getFlagSeverityColor(flag.severity)

    return (
        <div className={`rounded-lg border px-2.5 py-2 text-xs ${severityClasses}`}>
            <div className="flex items-center justify-between gap-2 cursor-pointer select-none" onClick={() => setOpen(!open)}>
                <div className="flex items-center gap-1.5 min-w-0">
                    {flag.severity === 'high' ? <AlertCircle className="w-3 h-3 flex-shrink-0" /> : <AlertTriangle className="w-3 h-3 flex-shrink-0" />}
                    <span className="font-medium truncate">{flag.property_name}</span>
                    <span className="font-mono opacity-70">{flag.value} {flag.direction} {flag.threshold}</span>
                </div>
                <span className={`flex-shrink-0 px-1.5 py-0.5 rounded text-[10px] font-semibold uppercase border ${severityClasses}`}>{flag.severity}</span>
            </div>
            {open && (
                <div className="mt-2 pt-2 border-t border-current/20 space-y-1.5 animate-slide-up">
                    <div className="flex items-start gap-1.5">
                        <AlertTriangle className="w-3 h-3 mt-0.5 flex-shrink-0 opacity-60" />
                        <p className="opacity-80">{flag.implication}</p>
                    </div>
                    <div className="flex items-start gap-1.5">
                        <Lightbulb className="w-3 h-3 mt-0.5 flex-shrink-0 opacity-60" />
                        <p className="opacity-70">{flag.recommendation}</p>
                    </div>
                </div>
            )}
        </div>
    )
}

// ── Download Button ───────────────────────────────────────────────────────────

interface DownloadBtnProps {
    label: string
    fileType: DockingFileType
    tooltip: string
    jobId: string
    compoundIndex: number
}

function DownloadBtn({ label, fileType, tooltip, jobId, compoundIndex }: DownloadBtnProps) {
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState<string | null>(null)
    const [tooltipOpen, setTooltipOpen] = useState(false)

    const handleClick = async () => {
        setLoading(true); setError(null)
        try { await downloadDockingFile(jobId, compoundIndex, fileType) }
        catch (e: unknown) {
            const msg = e instanceof Error ? e.message : 'Download failed'
            setError(msg)
            setTimeout(() => setError(null), 5000)
        } finally { setLoading(false) }
    }

    return (
        <div className="relative" onMouseEnter={() => !error && setTooltipOpen(true)} onMouseLeave={() => setTooltipOpen(false)}>
            <button
                onClick={handleClick}
                onFocus={() => !error && setTooltipOpen(true)}
                onBlur={() => setTooltipOpen(false)}
                disabled={loading}
                className={`flex items-center gap-1.5 px-2.5 py-1.5 rounded-md text-xs font-medium border transition-colors ${error
                    ? 'border-red-700 bg-red-950/40 text-red-400 hover:bg-red-900/40'
                    : 'border-gray-700 bg-gray-800/60 text-gray-400 hover:bg-gray-700 hover:text-gray-200'
                    } disabled:opacity-60 disabled:cursor-not-allowed`}
            >
                {loading ? <Loader2 className="w-3 h-3 animate-spin" /> : <Download className="w-3 h-3" />}
                <span>{label}</span>
            </button>
            {error
                ? <div className="absolute bottom-full left-0 mb-1.5 w-64 bg-gray-900 border border-red-800 rounded-lg px-2.5 py-2 text-[11px] text-red-400 z-50 shadow-xl leading-relaxed">{error}</div>
                : tooltipOpen
                    ? <span className="absolute bottom-full left-0 mb-1.5 w-60 bg-gray-900 border border-gray-700 rounded-lg px-2.5 py-2 text-xs text-gray-300 z-50 pointer-events-none leading-relaxed shadow-xl">{tooltip}</span>
                    : null}
        </div>
    )
}

// ── Narrative Block ───────────────────────────────────────────────────────────

function NarrativeCard({ block }: { block: ADMETNarrativeBlock }) {
    const style = getNarrativeStyle(block.level)
    return (
        <div className={`rounded-lg border px-3 py-2.5 ${style.wrap}`}>
            <div className="flex items-start gap-2">
                {getNarrativeIcon(block.level)}
                <div className="min-w-0">
                    <p className={`text-xs font-semibold mb-0.5 ${style.title}`}>{block.title}</p>
                    <p className={`text-xs leading-relaxed ${style.body}`}>{block.body}</p>
                </div>
            </div>
        </div>
    )
}

// ── Endpoint Table ────────────────────────────────────────────────────────────

function EndpointTable({ endpoints }: { endpoints: ADMETEndpointValue[] }) {
    const [activeCategory, setActiveCategory] = useState<string>('all')

    const { presentCategories, categoryCounts } = useMemo(() => {
        const counts: Record<string, number> = { all: endpoints.length }
        endpoints.forEach(e => { counts[e.category] = (counts[e.category] ?? 0) + 1 })
        const cats = new Set(endpoints.map(e => e.category))
        return { presentCategories: ['all', ...CATEGORY_ORDER.filter(c => cats.has(c))], categoryCounts: counts }
    }, [endpoints])

    const filtered = useMemo(
        () => activeCategory === 'all' ? endpoints : endpoints.filter(e => e.category === activeCategory),
        [endpoints, activeCategory],
    )

    return (
        <div>
            <div className="flex flex-wrap gap-1.5 mb-3">
                {presentCategories.map(cat => (
                    <button key={cat} onClick={() => setActiveCategory(cat)}
                        className={`px-2 py-0.5 rounded-full text-[10px] font-medium uppercase tracking-wider border transition-colors flex items-center gap-1 ${activeCategory === cat ? 'bg-gray-700 border-gray-500 text-gray-200' : 'border-gray-700 text-gray-500 hover:border-gray-600 hover:text-gray-400'}`}>
                        {cat === 'all' ? 'All' : CATEGORY_LABELS[cat]}
                        <span className={`text-[9px] ${activeCategory === cat ? 'text-gray-400' : 'text-gray-700'}`}>{categoryCounts[cat] ?? 0}</span>
                    </button>
                ))}
            </div>
            <div className="overflow-x-auto rounded-lg border border-gray-700">
                <table className="w-full text-xs">
                    <thead>
                        <tr className="border-b border-gray-700 bg-gray-800/80">
                            <th className="text-left px-3 py-2 text-gray-500 font-medium uppercase tracking-wider">Endpoint</th>
                            <th className="text-right px-3 py-2 text-gray-500 font-medium uppercase tracking-wider">Value</th>
                            <th className="text-left px-3 py-2 text-gray-500 font-medium uppercase tracking-wider hidden sm:table-cell">Threshold</th>
                            <th className="text-center px-3 py-2 text-gray-500 font-medium uppercase tracking-wider">Status</th>
                            <th className="text-left px-3 py-2 text-gray-500 font-medium uppercase tracking-wider hidden md:table-cell">Category</th>
                        </tr>
                    </thead>
                    <tbody>
                        {filtered.map(ep => (
                            <tr key={ep.key} className={`border-b border-gray-800 last:border-0 transition-colors ${ep.triggered ? ep.severity === 'high' ? 'bg-red-950/20 hover:bg-red-950/30' : 'bg-yellow-950/10 hover:bg-yellow-950/20' : 'hover:bg-gray-800/40'}`} title={ep.interpretation ?? ''}>
                                <td className="px-3 py-2">
                                    <div className="flex items-center gap-1.5">
                                        <span className="text-gray-300">{ep.label}</span>
                                        {ep.hard_fail && <span className="text-[9px] uppercase tracking-wider text-red-500 border border-red-900 rounded px-1 py-0.5">hard fail</span>}
                                        {ep.interpretation && <TooltipBubble text={ep.interpretation} widthClass="w-60" align="left" iconClassName="w-3 h-3 text-gray-700 hover:text-gray-500" className="ml-0" />}
                                    </div>
                                </td>
                                <td className="px-3 py-2 text-right font-mono">
                                    <span className={getEndpointStatusClass(ep)}>{ep.display_value}</span>
                                    {ep.unit && <span className="text-gray-600 ml-1 text-[10px]">{ep.unit}</span>}
                                </td>
                                <td className="px-3 py-2 text-gray-600 font-mono hidden sm:table-cell">{ep.threshold_applied ?? '—'}</td>
                                <td className="px-3 py-2 text-center">
                                    {ep.threshold_applied
                                        ? ep.triggered
                                            ? <span className={`inline-flex items-center gap-1 text-[10px] font-medium ${ep.severity === 'high' ? 'text-red-400' : 'text-yellow-400'}`}><AlertTriangle className="w-3 h-3" />{ep.severity}</span>
                                            : <span className="inline-flex items-center gap-1 text-[10px] font-medium text-emerald-400"><CheckCircle2 className="w-3 h-3" />pass</span>
                                        : <span className="text-gray-600 text-[10px]">—</span>}
                                </td>
                                <td className="px-3 py-2 hidden md:table-cell">
                                    <span className="text-[10px] uppercase tracking-wider text-gray-600 border border-gray-800 rounded px-1.5 py-0.5">{CATEGORY_LABELS[ep.category] ?? ep.category}</span>
                                </td>
                            </tr>
                        ))}
                    </tbody>
                </table>
            </div>
        </div>
    )
}

// ── Active Thresholds Panel ───────────────────────────────────────────────────

function ActiveThresholdsPanel({ thresholds }: { thresholds: Record<string, ADMETThresholdConfig> }) {
    const entries = Object.entries(thresholds)
    if (entries.length === 0) return null
    return (
        <div className="overflow-x-auto rounded-lg border border-gray-700/60">
            <table className="w-full text-xs">
                <thead>
                    <tr className="border-b border-gray-700/60 bg-gray-800/40">
                        <th className="text-left px-3 py-1.5 text-gray-600 font-medium uppercase tracking-wider">Endpoint</th>
                        <th className="text-left px-3 py-1.5 text-gray-600 font-medium uppercase tracking-wider">Cutoff</th>
                        <th className="text-center px-3 py-1.5 text-gray-600 font-medium uppercase tracking-wider">Hard Fail</th>
                        <th className="text-center px-3 py-1.5 text-gray-600 font-medium uppercase tracking-wider">Active</th>
                    </tr>
                </thead>
                <tbody>
                    {entries.map(([key, cfg]) => (
                        <tr key={key} className="border-b border-gray-800/60 last:border-0 hover:bg-gray-800/30">
                            <td className="px-3 py-1.5 text-gray-400">{THRESHOLD_KEY_LABELS[key] ?? key}</td>
                            <td className="px-3 py-1.5 font-mono text-gray-400">{cfg.threshold_str ?? `${cfg.direction === 'above' ? '>' : '<'}${cfg.cutoff}`}</td>
                            <td className="px-3 py-1.5 text-center">
                                {cfg.hard_fail ? <span className="text-[10px] text-red-500 border border-red-900/60 rounded px-1.5 py-0.5 uppercase tracking-wider">yes</span> : <span className="text-gray-700 text-[10px]">—</span>}
                            </td>
                            <td className="px-3 py-1.5 text-center">
                                {cfg.enabled ? <CheckCircle2 className="w-3 h-3 text-emerald-600 inline" /> : <span className="text-gray-700 text-[10px]">off</span>}
                            </td>
                        </tr>
                    ))}
                </tbody>
            </table>
        </div>
    )
}

// ── Full ADMET Report ─────────────────────────────────────────────────────────

function FullADMETReport({ admet }: { admet: ExtendedADMET }) {
    const [endpointTableOpen, setEndpointTableOpen] = useState(false)
    const [thresholdsOpen, setThresholdsOpen] = useState(false)

    const hasNarrative = (admet.narrative?.length ?? 0) > 0
    const hasEndpointTable = (admet.endpoint_table?.length ?? 0) > 0
    const hasDecisionBasis = (admet.decision_basis?.length ?? 0) > 0
    const hasThresholds = Object.keys(admet.thresholds_used ?? {}).length > 0

    if (!hasNarrative && !hasEndpointTable && !hasDecisionBasis) return null

    return (
        <div className="rounded-lg bg-gray-800/30 border border-gray-700 overflow-hidden">
            <div className="flex items-center gap-2 px-3 py-2.5 border-b border-gray-700 bg-gray-800/50">
                <BookOpen className="w-3.5 h-3.5 text-gray-500" />
                <p className="text-xs font-medium text-gray-300 uppercase tracking-wider">Full ADMET Report</p>
                <span className="text-[10px] text-gray-600 ml-auto">ADMET-AI · 41 TDC datasets</span>
            </div>
            <div className="p-3 space-y-4">
                {hasNarrative && <div className="space-y-2">{admet.narrative!.map((block, i) => <NarrativeCard key={i} block={block} />)}</div>}
                {hasEndpointTable && (
                    <div>
                        <button onClick={() => setEndpointTableOpen(v => !v)} className="w-full flex items-center justify-between gap-2 text-xs text-gray-400 hover:text-gray-300 transition-colors py-1">
                            <div className="flex items-center gap-1.5">
                                <Layers className="w-3.5 h-3.5 text-gray-500" />
                                <span className="font-medium uppercase tracking-wider">Full Endpoint Table</span>
                                <span className="text-gray-600">{admet.endpoint_table!.length} properties</span>
                            </div>
                            {endpointTableOpen ? <ChevronUp className="w-3.5 h-3.5 text-gray-600" /> : <ChevronDown className="w-3.5 h-3.5 text-gray-600" />}
                        </button>
                        {endpointTableOpen && <div className="mt-2 animate-slide-up"><EndpointTable endpoints={admet.endpoint_table!} /></div>}
                    </div>
                )}
                {hasThresholds && (
                    <div>
                        <button onClick={() => setThresholdsOpen(v => !v)} className="w-full flex items-center justify-between gap-2 text-xs text-gray-400 hover:text-gray-300 transition-colors py-1">
                            <div className="flex items-center gap-1.5">
                                <SlidersHorizontal className="w-3.5 h-3.5 text-gray-500" />
                                <span className="font-medium uppercase tracking-wider">Active Thresholds</span>
                                <span className="text-gray-600 text-[10px]">applied to this result</span>
                            </div>
                            {thresholdsOpen ? <ChevronUp className="w-3.5 h-3.5 text-gray-600" /> : <ChevronDown className="w-3.5 h-3.5 text-gray-600" />}
                        </button>
                        {thresholdsOpen && <div className="mt-2 animate-slide-up"><ActiveThresholdsPanel thresholds={admet.thresholds_used!} /></div>}
                    </div>
                )}
                {hasDecisionBasis && (
                    <div>
                        <p className="text-[10px] text-gray-600 uppercase tracking-wider mb-1.5">Decision Basis</p>
                        <ul className="space-y-0.5">
                            {admet.decision_basis!.map((reason, i) => (
                                <li key={i} className="flex items-start gap-1.5 text-xs text-gray-600">
                                    <span className="mt-1.5 w-1 h-1 rounded-full bg-gray-700 flex-shrink-0" />
                                    {reason}
                                </li>
                            ))}
                        </ul>
                    </div>
                )}
            </div>
        </div>
    )
}

// ── Retrosynthesis Step Card ──────────────────────────────────────────────────

function RetroStepCard({ step }: { step: ExtendedRetrosynthesisStep }) {
    const [open, setOpen] = useState(false)
    const [protocolOpen, setProtocolOpen] = useState(false)

    const hasDetails = !!(step.reaction_name || step.reagents?.length || step.conditions || step.temperature || step.duration || step.yield_estimate || step.protocol_text)

    return (
        <div className="rounded-lg border border-gray-700 bg-gray-900/40 overflow-hidden">
            <div className={`flex items-center justify-between gap-3 px-3 py-2.5 ${hasDetails ? 'cursor-pointer hover:bg-gray-800/40 transition-colors' : ''}`} onClick={() => hasDetails && setOpen(!open)}>
                <div className="flex items-center gap-2 min-w-0">
                    <span className="flex-shrink-0 w-5 h-5 rounded-full bg-gray-700 text-gray-300 text-[10px] font-bold flex items-center justify-center">{step.step_number}</span>
                    <div className="min-w-0">
                        {step.reaction_name ? <p className="text-xs font-medium text-gray-200 truncate">{step.reaction_name}</p> : <p className="text-xs text-gray-500 italic">Unnamed transformation</p>}
                        {step.reaction_type && <p className="text-[10px] text-gray-600 mt-0.5">{REACTION_TYPE_LABELS[step.reaction_type] ?? step.reaction_type}</p>}
                    </div>
                </div>
                <div className="flex items-center gap-2 flex-shrink-0">
                    {step.yield_estimate && (
                        <span className="flex items-center gap-0.5 text-[10px] text-emerald-500 bg-emerald-950/40 border border-emerald-900/60 rounded px-1.5 py-0.5">
                            <Percent className="w-2.5 h-2.5" />{step.yield_estimate}
                        </span>
                    )}
                    <span className="text-[10px] text-gray-600">{(step.confidence * 100).toFixed(0)}% conf.</span>
                    {hasDetails && (open ? <ChevronUp className="w-3.5 h-3.5 text-gray-600" /> : <ChevronDown className="w-3.5 h-3.5 text-gray-600" />)}
                </div>
            </div>

            {step.schematic && (
                <div className="px-3 pb-2.5 border-t border-gray-800/60">
                    <div className="flex items-center gap-1.5 mt-2 overflow-x-auto">
                        {step.schematic.split('→').map((part, i, arr) => (
                            <div key={i} className="flex items-center gap-1.5 flex-shrink-0">
                                <code className="font-mono text-[10px] text-teal-400 bg-gray-800 px-2 py-1 rounded border border-gray-700 max-w-[180px] truncate block" title={part.trim()}>{part.trim()}</code>
                                {i < arr.length - 1 && <ArrowRight className="w-3 h-3 text-gray-600 flex-shrink-0" />}
                            </div>
                        ))}
                    </div>
                </div>
            )}

            {open && hasDetails && (
                <div className="px-3 pb-3 border-t border-gray-800 space-y-3 pt-3 animate-slide-up">
                    {(step.temperature || step.duration || step.conditions) && (
                        <div className="grid grid-cols-2 gap-2 sm:grid-cols-3">
                            {step.temperature && (
                                <div className="flex items-center gap-1.5 text-xs">
                                    <Thermometer className="w-3 h-3 text-gray-600 flex-shrink-0" />
                                    <div><p className="text-[10px] text-gray-600 uppercase tracking-wider">Temp</p><p className="text-gray-300">{step.temperature}</p></div>
                                </div>
                            )}
                            {step.duration && (
                                <div className="flex items-center gap-1.5 text-xs">
                                    <Clock className="w-3 h-3 text-gray-600 flex-shrink-0" />
                                    <div><p className="text-[10px] text-gray-600 uppercase tracking-wider">Duration</p><p className="text-gray-300">{step.duration}</p></div>
                                </div>
                            )}
                            {step.yield_estimate && (
                                <div className="flex items-center gap-1.5 text-xs">
                                    <Percent className="w-3 h-3 text-gray-600 flex-shrink-0" />
                                    <div><p className="text-[10px] text-gray-600 uppercase tracking-wider">Est. Yield</p><p className="text-emerald-400">{step.yield_estimate}</p></div>
                                </div>
                            )}
                        </div>
                    )}
                    {step.conditions && <div><p className="text-[10px] text-gray-600 uppercase tracking-wider mb-1">Conditions</p><p className="text-xs text-gray-400 leading-relaxed">{step.conditions}</p></div>}
                    {(step.reagents?.length ?? 0) > 0 && (
                        <div>
                            <div className="flex items-center gap-1.5 mb-1.5"><FlaskConical className="w-3 h-3 text-gray-600" /><p className="text-[10px] text-gray-600 uppercase tracking-wider">Reagents</p></div>
                            <div className="flex flex-wrap gap-1">{step.reagents!.map((r, i) => <span key={i} className="text-[10px] text-gray-400 bg-gray-800 border border-gray-700 rounded px-1.5 py-0.5 font-mono">{r}</span>)}</div>
                        </div>
                    )}
                    {(step.solvents?.length ?? 0) > 0 && (
                        <div>
                            <div className="flex items-center gap-1.5 mb-1.5"><TestTube className="w-3 h-3 text-gray-600" /><p className="text-[10px] text-gray-600 uppercase tracking-wider">Solvents</p></div>
                            <div className="flex flex-wrap gap-1">{step.solvents!.map((s, i) => <span key={i} className="text-[10px] text-blue-400/80 bg-blue-950/30 border border-blue-900/50 rounded px-1.5 py-0.5 font-mono">{s}</span>)}</div>
                        </div>
                    )}
                    {step.starting_materials.length > 0 && (
                        <div>
                            <p className="text-[10px] text-gray-600 uppercase tracking-wider mb-1.5">Starting Materials</p>
                            <div className="space-y-1">{step.starting_materials.map((smi, i) => <code key={i} className="block font-mono text-[10px] text-teal-400 bg-gray-800 px-2 py-1 rounded border border-gray-700 truncate" title={smi}>{smi}</code>)}</div>
                        </div>
                    )}
                    {step.protocol_text && (
                        <div className="border border-gray-700 rounded-lg overflow-hidden">
                            <button onClick={e => { e.stopPropagation(); setProtocolOpen(!protocolOpen) }} className="w-full flex items-center justify-between gap-2 px-3 py-2 text-xs text-gray-400 hover:text-gray-300 hover:bg-gray-800/40 transition-colors">
                                <div className="flex items-center gap-1.5"><FileText className="w-3.5 h-3.5 text-gray-600" /><span className="font-medium uppercase tracking-wider text-[10px]">Lab Protocol</span></div>
                                {protocolOpen ? <ChevronUp className="w-3.5 h-3.5 text-gray-600" /> : <ChevronDown className="w-3.5 h-3.5 text-gray-600" />}
                            </button>
                            {protocolOpen && <div className="px-3 pb-3 pt-1 border-t border-gray-700 animate-slide-up"><p className="text-xs text-gray-400 leading-relaxed whitespace-pre-line">{step.protocol_text}</p></div>}
                        </div>
                    )}
                </div>
            )}
        </div>
    )
}

// ── Retrosynthesis Panel ──────────────────────────────────────────────────────

function RetrosynthesisPanel({ retro }: { retro: ExtendedRetrosynthesisResult }) {
    const difficultyColor = DIFFICULTY_STYLES[retro.difficulty_label ?? 'unknown'] ?? 'text-gray-500'
    const complexityInfo = getComplexityLabel(retro.complexity_score)

    return (
        <div className="p-3 rounded-lg bg-gray-800/50 border border-gray-700 overflow-visible">
            <div className="flex items-center gap-1.5 mb-3">
                <GitBranch className="w-3.5 h-3.5 text-gray-500" />
                <p className="text-xs font-medium text-gray-300 uppercase tracking-wider">Retrosynthesis</p>
                <Tip text="Analyses whether this molecule can be synthesized in a lab using BRICS retrosynthetic fragmentation. SA Score rates synthetic accessibility 1 (trivial) → 10 (infeasible)." />
            </div>
            {retro.feasible ? (
                <div className="space-y-3">
                    <div className="grid grid-cols-2 gap-x-4 gap-y-2">
                        <div className="flex justify-between text-xs"><span className="text-gray-500">Synthesis steps</span><span className="text-gray-200 font-medium">{retro.num_steps}</span></div>
                        {retro.sa_score !== undefined && retro.sa_score !== null && retro.sa_score > 0 && (
                            <div className="flex justify-between text-xs">
                                <span className="text-gray-500 flex items-center">SA Score <Tip text="1.0 = trivially easy, 10.0 = practically impossible. ≤3.5: easy, 3.5–6.0: feasible, >6.0: hard." /></span>
                                <span className={difficultyColor}>{retro.sa_score.toFixed(2)}</span>
                            </div>
                        )}
                        <div className="flex justify-between text-xs items-center">
                            <span className="text-gray-500 flex items-center">Difficulty <Tip text="Easy: ≤3.5, Moderate: ≤5.0, Hard: ≤6.0, Infeasible: >6.0." /></span>
                            <span className={`font-medium capitalize ${difficultyColor}`}>{retro.difficulty_label ?? getComplexityLabel(retro.complexity_score).label}</span>
                        </div>
                        <div className="flex justify-between text-xs">
                            <span className="text-gray-500 flex items-center">Complexity <Tip text="SA Score × 10. ≤15 Very Easy, 15–25 Moderate, 25–40 Difficult, >40 Infeasible." /></span>
                            <span className={complexityInfo.color}>{retro.complexity_score.toFixed(1)}</span>
                        </div>
                        {retro.estimated_total_yield && (
                            <div className="flex justify-between text-xs col-span-2">
                                <span className="text-gray-500 flex items-center">Estimated total yield <Tip text="Multiplicative estimate across all steps. Rough approximation only." /></span>
                                <span className="text-emerald-400 font-medium">{retro.estimated_total_yield}</span>
                            </div>
                        )}
                    </div>
                    {retro.synthesis_summary && (
                        <div className="rounded-lg bg-gray-800/60 border border-gray-700 px-3 py-2.5">
                            <p className="text-[10px] text-gray-600 uppercase tracking-wider mb-1">Summary</p>
                            <p className="text-xs text-gray-400 leading-relaxed">{retro.synthesis_summary}</p>
                        </div>
                    )}
                    {retro.route.length > 0 && (
                        <div className="space-y-2">
                            <p className="text-[10px] text-gray-600 uppercase tracking-wider">Synthesis Route · {retro.route.length} step{retro.route.length !== 1 ? 's' : ''}</p>
                            {retro.route.map(step => <RetroStepCard key={step.step_number} step={step} />)}
                        </div>
                    )}
                </div>
            ) : (
                <div className="rounded-lg bg-red-950/20 border border-red-900/50 px-3 py-2.5">
                    <div className="flex items-start gap-2">
                        <AlertCircle className="w-3.5 h-3.5 text-red-500 flex-shrink-0 mt-0.5" />
                        <div>
                            <p className="text-xs font-medium text-red-400 mb-0.5">No feasible synthesis route</p>
                            {retro.synthesis_summary ? <p className="text-xs text-red-300/60 leading-relaxed">{retro.synthesis_summary}</p> : <p className="text-xs text-red-300/60">SA Score too high — molecular complexity exceeds practical synthesis limits.</p>}
                            {retro.sa_score !== undefined && retro.sa_score !== null && retro.sa_score > 0 && <p className="text-[10px] text-red-500/70 mt-1 font-mono">SA Score: {retro.sa_score.toFixed(2)}</p>}
                        </div>
                    </div>
                </div>
            )}
        </div>
    )
}

// ── Score Breakdown Parser ────────────────────────────────────────────────────
// Reads the v4 breakdown dict from backend.
// Scored rows: have contribution + max_possible (= effective_weight, a whole integer)
// Gate-only row: binding_prefilter with gate_only:true — rendered as info pill
// Info-only rows: weight_redistribution, penalties, final_score — skipped from scoring display

interface ParsedBreakdown {
    scoredRows: Array<{
        key: string
        entry: BreakdownEntry
        contribution: number   // pts earned (float, e.g. 25.5)
        maxPossible: number    // whole-number max pts (effective_weight)
        raw: string
        pct: number            // 0-100
    }>
    gateRow: BreakdownEntry | null    // binding_prefilter info
    totalMax: number                   // sum of all active maxPossible
    finalScore: number | null
    penaltyApplied: boolean
    disabledComponents: string[]
}

function parseBreakdown(bd: ScoreBreakdown): ParsedBreakdown {
    const raw = bd as Record<string, unknown>

    const scoredRows: ParsedBreakdown['scoredRows'] = []
    let gateRow: BreakdownEntry | null = null
    let finalScore: number | null = null
    let penaltyApplied = false
    let disabledComponents: string[] = []

    // final_score
    if (typeof raw['final_score'] === 'number') {
        finalScore = raw['final_score'] as number
    }

    // penalties
    const penalties = raw['penalties'] as Record<string, unknown> | undefined
    if (penalties) {
        penaltyApplied = Boolean(penalties['mw_fragment_penalty_applied'])
    }

    // weight_redistribution
    const wr = raw['weight_redistribution'] as Record<string, unknown> | undefined
    if (wr?.disabled_components && Array.isArray(wr.disabled_components)) {
        disabledComponents = wr.disabled_components as string[]
    }

    // binding_prefilter gate row
    const bpRaw = raw['binding_prefilter'] as BreakdownEntry | undefined
    if (bpRaw?.gate_only) {
        gateRow = bpRaw
    }

    // Scored rows — iterate in defined order, then any extra keys
    const allKeys = Object.keys(raw)
    const orderedKeys = [
        ...SCORED_KEY_ORDER.filter(k => allKeys.includes(k)),
        ...allKeys.filter(k => !SCORED_KEY_ORDER.includes(k) && !INFO_ONLY_KEYS.has(k) && k !== 'binding_prefilter'),
    ]

    for (const key of orderedKeys) {
        const entry = raw[key] as BreakdownEntry | undefined
        if (!entry || typeof entry !== 'object') continue
        if (entry.gate_only) continue                         // prefilter handled above
        if (INFO_ONLY_KEYS.has(key)) continue

        const contribution = entry.contribution ?? 0
        // max_possible = effective_weight (whole-number pts from backend v4)
        const maxPossible = entry.max_possible ?? entry.effective_weight ?? 0
        const rawStr = (entry.raw ?? '').trim()
        const pct = maxPossible > 0 ? Math.min((contribution / maxPossible) * 100, 100) : 0

        // Skip rows where step is disabled (effective_weight = 0) and no contribution
        if (maxPossible === 0 && contribution === 0 && entry.step_active === false) continue

        scoredRows.push({ key, entry, contribution, maxPossible, raw: rawStr, pct })
    }

    const totalMax = scoredRows
        .filter(r => r.maxPossible > 0)
        .reduce((s, r) => s + r.maxPossible, 0)

    return { scoredRows, gateRow, totalMax, finalScore, penaltyApplied, disabledComponents }
}

// ── Main Card ─────────────────────────────────────────────────────────────────

interface MoleculeCardProps {
    compound: CompoundResult
    jobId: string
    index: number
}

export default function MoleculeCard({ compound, jobId, index }: MoleculeCardProps) {
    const [expanded, setExpanded] = useState(false)
    const [copied, setCopied] = useState(false)
    const [breakdown, setBreakdown] = useState<ScoreBreakdown | null>(null)
    const [breakdownLoading, setBreakdownLoading] = useState(false)

    const score = compound.final_score ?? 0
    const admet = compound.admet as ExtendedADMET | null | undefined
    const retro = compound.retrosynthesis as ExtendedRetrosynthesisResult | null | undefined

    const copySmiles = useCallback(() => {
        try {
            if (navigator?.clipboard?.writeText) {
                navigator.clipboard.writeText(compound.canonical_smiles).catch(() => {
                    const el = document.createElement('textarea')
                    el.value = compound.canonical_smiles
                    el.style.cssText = 'position:fixed;opacity:0'
                    document.body.appendChild(el); el.select(); document.execCommand('copy'); document.body.removeChild(el)
                })
            }
            setCopied(true); setTimeout(() => setCopied(false), 2000)
        } catch { }
    }, [compound.canonical_smiles])

    const handleExpand = useCallback(() => {
        const willExpand = !expanded
        setExpanded(willExpand)
        if (willExpand && !breakdown) {
            setBreakdownLoading(true)
            getScoreBreakdown(jobId, index)
                .then(data => setBreakdown(data))
                .catch(() => { })
                .finally(() => setBreakdownLoading(false))
        }
    }, [expanded, breakdown, jobId, index])

    const ringStyle = score >= 70 ? 'border-emerald-500 text-emerald-400' : score >= 45 ? 'border-yellow-500 text-yellow-400' : 'border-red-600 text-red-400'
    const flagCount = admet?.flags?.length ?? admet?.flag_summary?.length ?? 0

    const parsedBreakdown = useMemo(
        () => breakdown ? parseBreakdown(breakdown) : null,
        [breakdown],
    )

    return (
        <div className="card-hover animate-fade-in overflow-visible">
            {/* ── Summary row ── */}
            <div className="flex items-center gap-4">
                <div className="flex-shrink-0 text-center">
                    <div className={`score-ring ${ringStyle}`}>{score.toFixed(0)}</div>
                    <p className="text-xs text-gray-600 mt-1">#{compound.rank}</p>
                </div>
                <div className="flex-1 min-w-0">
                    <div className="flex items-center gap-2">
                        <span className="smiles-display">{compound.canonical_smiles}</span>
                        <button onClick={copySmiles} className="text-gray-600 hover:text-gray-400 transition-colors flex-shrink-0" title="Copy SMILES">
                            {copied ? <CheckCheck className="w-3.5 h-3.5 text-emerald-500" /> : <Copy className="w-3.5 h-3.5" />}
                        </button>
                    </div>
                    <div className="flex flex-wrap items-center gap-x-4 gap-y-1 mt-2">
                        {compound.lipinski && <span className={compound.lipinski.passed ? 'badge-pass' : 'badge-fail'}>Lipinski {compound.lipinski.passed ? '✓' : '✗'}</span>}
                        {admet && <span className={admet.passed ? 'badge-pass' : 'badge-warn'}>ADMET {flagCount > 0 ? `${flagCount} flag${flagCount !== 1 ? 's' : ''}` : 'clean'}</span>}
                        {admet?.overall_risk && admet.overall_risk !== 'low' && (
                            <span className={`inline-flex items-center gap-1 text-[10px] font-semibold uppercase px-1.5 py-0.5 rounded border ${getRiskStyle(admet.overall_risk)}`}>
                                <Shield className="w-2.5 h-2.5" />{admet.overall_risk} risk
                            </span>
                        )}
                        {compound.docking && (
                            <span className={`text-xs font-mono font-medium ${getAffinityColor(compound.docking.best_affinity_kcal)}`}>
                                {compound.docking.best_affinity_kcal.toFixed(2)} kcal/mol
                            </span>
                        )}
                        {retro && (
                            <span className={`text-xs ${retro.feasible ? 'text-gray-500' : 'text-red-500/70'}`}>
                                {retro.feasible ? `${retro.num_steps} synthesis step${retro.num_steps !== 1 ? 's' : ''}${retro.difficulty_label ? ` · ${retro.difficulty_label}` : ''}` : 'synthesis infeasible'}
                            </span>
                        )}
                    </div>
                </div>
                <button onClick={handleExpand} className="flex-shrink-0 text-gray-600 hover:text-gray-400 transition-colors p-1" aria-label={expanded ? 'Collapse' : 'Expand'}>
                    {expanded ? <ChevronUp className="w-4 h-4" /> : <ChevronDown className="w-4 h-4" />}
                </button>
            </div>

            {/* ── Expanded detail ── */}
            {expanded && (
                <div className="mt-4 pt-4 border-t border-gray-800 space-y-5 animate-slide-up overflow-visible">

                    {/* SMILES */}
                    <div>
                        <p className="text-xs text-gray-500 mb-1 uppercase tracking-wider">Canonical SMILES</p>
                        <p className="font-mono text-xs text-emerald-400 bg-gray-800 px-3 py-2 rounded-lg border border-gray-700 break-all">{compound.canonical_smiles}</p>
                    </div>

                    {/* Drug-likeness + ADMET side by side */}
                    <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
                        {compound.lipinski && (
                            <div className="p-3 rounded-lg bg-gray-800/50 border border-gray-700 overflow-visible">
                                <div className="flex items-center gap-1.5 mb-2">
                                    <FlaskConical className="w-3.5 h-3.5 text-gray-500" />
                                    <p className="text-xs font-medium text-gray-300 uppercase tracking-wider">Drug-likeness</p>
                                    <Tip text="Lipinski Rule of 5 — checks if this molecule has physical properties typical of orally absorbed drugs." />
                                </div>
                                <div className="space-y-1.5">
                                    {([
                                        { label: 'Mol. Weight', value: `${compound.lipinski.mw.toFixed(1)} Da`, ok: compound.lipinski.mw <= 500, tip: 'Ideal ≤500 Da. Heavier molecules struggle to cross cell membranes.' },
                                        { label: 'LogP', value: compound.lipinski.logp.toFixed(2), ok: compound.lipinski.logp <= 4.5, tip: 'Ideal 1–3 (optimal), 0–5 (ok). Measures fat-solubility. Scoring penalises >5.' },
                                        { label: 'H-bond donors', value: compound.lipinski.hbd, ok: compound.lipinski.hbd <= 5, tip: 'Ideal ≤5. -OH and -NH groups. Too many prevent gut absorption.' },
                                        { label: 'H-bond acceptors', value: compound.lipinski.hba, ok: compound.lipinski.hba <= 10, tip: 'Ideal ≤10. N and O atoms. Too many reduce oral bioavailability.' },
                                        { label: 'LogS solubility', value: compound.lipinski.logs.toFixed(2), ok: compound.lipinski.logs >= -4, tip: 'v4 scoring: ≥−3 Excellent (0.90), ≥−4 Good (0.75), ≥−5 Borderline (0.45), ≥−6 Poor (0.20), <−6 Insoluble (0.05).' },
                                    ] as const).map(({ label, value, ok, tip }) => (
                                        <div key={label} className="flex justify-between text-xs items-center">
                                            <span className="text-gray-500 flex items-center">{label}<Tip text={tip} /></span>
                                            <span className={ok ? 'text-emerald-400' : 'text-yellow-400'}>{value}</span>
                                        </div>
                                    ))}
                                </div>
                                <div className="flex justify-between text-xs pt-1 border-t border-gray-700 items-center mt-1">
                                    <span className="text-gray-500 flex items-center">Solubility class <Tip text="High: logS ≥ −3, Low: −5 to −3, Very Low: < −5. Class is informational; actual scoring uses continuous logS breakpoints." /></span>
                                    <span className={compound.lipinski.solubility_class === 'Very Low' ? 'text-red-400' : compound.lipinski.solubility_class === 'Low' ? 'text-yellow-400' : 'text-gray-300'}>{compound.lipinski.solubility_class}</span>
                                </div>
                            </div>
                        )}

                        {admet && (
                            <div className="p-3 rounded-lg bg-gray-800/50 border border-gray-700 overflow-visible">
                                <div className="flex items-center justify-between mb-2">
                                    <div className="flex items-center gap-1.5">
                                        <Activity className="w-3.5 h-3.5 text-gray-500" />
                                        <p className="text-xs font-medium text-gray-300 uppercase tracking-wider">ADMET Profile</p>
                                        <Tip text="Absorption, Distribution, Metabolism, Excretion, Toxicity — main reasons drugs fail clinical trials." />
                                    </div>
                                    {admet.overall_risk && (
                                        <span className={`inline-flex items-center gap-1 text-[10px] font-semibold uppercase px-1.5 py-0.5 rounded border ${getRiskStyle(admet.overall_risk)}`}>
                                            <Shield className="w-2.5 h-2.5" />{admet.overall_risk}
                                        </span>
                                    )}
                                </div>
                                <div className="space-y-1.5">
                                    {([
                                        { label: 'hERG inhibition', value: formatProbability(admet.herg_inhibition), risk: admet.herg_inhibition > 0.5, tip: 'Ideal <30%. Blocking hERG causes fatal arrhythmia — the #1 reason drugs are withdrawn post-market.' },
                                        { label: 'Hepatotoxicity', value: formatProbability(admet.hepatotoxicity), risk: admet.hepatotoxicity > 0.5, tip: 'Ideal <30%. Liver damage probability. #1 reason drugs fail FDA approval after Phase II.' },
                                        { label: 'Oral bioavailability', value: formatProbability(admet.oral_bioavailability), risk: admet.oral_bioavailability < 0.3, tip: 'Ideal >70%. Fraction of swallowed dose reaching bloodstream.' },
                                        { label: 'BBB penetration', value: formatProbability(admet.bbb_penetration), risk: false, tip: 'Desirable for CNS drugs; liability for non-CNS.' },
                                        { label: 'Caco-2 permeability', value: admet.caco2_permeability.toFixed(2), risk: admet.caco2_permeability < -5.15, tip: 'Ideal > −5.15. Models intestinal wall crossing. Key predictor of oral absorption.' },
                                    ] as const).map(({ label, value, risk, tip }) => (
                                        <div key={label} className="flex justify-between text-xs items-center">
                                            <span className="text-gray-500 flex items-center">{label}<Tip text={tip} /></span>
                                            <span className={risk ? 'text-red-400' : 'text-emerald-400'}>{value}</span>
                                        </div>
                                    ))}
                                </div>
                                {admet.flags.length > 0 && (
                                    <div className="mt-2 pt-2 border-t border-gray-700 space-y-1.5">
                                        <p className="text-[10px] text-gray-600 uppercase tracking-wider mb-1">Flags · click to expand</p>
                                        {admet.flags.map((flag, i) => <ADMETFlagCard key={i} flag={flag} />)}
                                    </div>
                                )}
                            </div>
                        )}
                    </div>

                    {/* Docking */}
                    {compound.docking && (
                        <div className="p-3 rounded-lg bg-gray-800/50 border border-gray-700 overflow-visible">
                            <div className="flex items-center gap-1.5 mb-2">
                                <Dna className="w-3.5 h-3.5 text-gray-500" />
                                <p className="text-xs font-medium text-gray-300 uppercase tracking-wider">Vina Docking</p>
                                <Tip text="AutoDock Vina simulates the drug fitting into the target protein's active site. More negative = stronger binding." />
                            </div>
                            <div className="space-y-1.5">
                                <div className="flex justify-between text-xs items-center">
                                    <span className="text-gray-500 flex items-center">Best affinity <Tip text="≤−9.0: Outstanding, −7 to −8.9: Strong, −5 to −6.9: Moderate, >−5: Weak." /></span>
                                    <div className="flex items-center gap-2">
                                        <span className={`font-mono font-medium ${getAffinityColor(compound.docking.best_affinity_kcal)}`}>{compound.docking.best_affinity_kcal.toFixed(3)} kcal/mol</span>
                                        <span className={`text-xs ${getAffinityLabel(compound.docking.best_affinity_kcal).color}`}>{getAffinityLabel(compound.docking.best_affinity_kcal).label}</span>
                                    </div>
                                </div>
                                <div className="flex justify-between text-xs"><span className="text-gray-500">Poses generated</span><span className="text-gray-300">{compound.docking.poses.length}</span></div>
                            </div>
                            {compound.docking.poses.length > 1 && (
                                <div className="mt-2 pt-2 border-t border-gray-700">
                                    <p className="text-xs text-gray-600 mb-1">Top poses</p>
                                    <div className="space-y-0.5">
                                        {compound.docking.poses.slice(0, 4).map(pose => (
                                            <div key={pose.rank} className="flex justify-between text-xs font-mono">
                                                <span className="text-gray-600">#{pose.rank}</span>
                                                <span className="text-gray-400">{pose.affinity_kcal.toFixed(2)}</span>
                                                <span className="text-gray-600">rmsd {pose.rmsd_ub.toFixed(2)}</span>
                                            </div>
                                        ))}
                                    </div>
                                </div>
                            )}
                            <div className="mt-3 pt-3 border-t border-gray-700">
                                <p className="text-[10px] text-gray-600 uppercase tracking-wider mb-2 flex items-center gap-1"><Download className="w-3 h-3" /> Download docking files</p>
                                <div className="flex flex-wrap gap-2">
                                    <DownloadBtn label="Raw .pdbqt" fileType="pdbqt" jobId={jobId} compoundIndex={index} tooltip="All Vina poses. Load in PyMOL or pass to PLIP manually." />
                                    <DownloadBtn label="Pose 1+H .pdb" fileType="pose_pdb" jobId={jobId} compoundIndex={index} tooltip="Best pose only, hydrogens added via obabel. Ready for PLIP or PyMOL." />
                                    <DownloadBtn label="PLIP Complex .pdb" fileType="complex_pdb" jobId={jobId} compoundIndex={index} tooltip="Receptor + docked ligand merged. Drop directly into PLIP — no manual merging." />
                                </div>
                            </div>
                        </div>
                    )}

                    {retro && <RetrosynthesisPanel retro={retro} />}
                    {admet && <FullADMETReport admet={admet} />}

                    {/* Score Breakdown */}
                    {breakdownLoading && (
                        <div className="flex items-center gap-2 text-xs text-gray-600">
                            <div className="w-3 h-3 border border-gray-600 border-t-gray-400 rounded-full animate-spin" />
                            Loading score breakdown...
                        </div>
                    )}

                    {parsedBreakdown && !breakdownLoading && (
                        <div className="overflow-visible">
                            {/* Header */}
                            <div className="flex items-center justify-between mb-3">
                                <p className="text-xs text-gray-500 uppercase tracking-wider">Score Breakdown</p>
                                <span className="text-xs font-mono">
                                    <span className={score >= 70 ? 'text-emerald-400' : score >= 45 ? 'text-yellow-400' : 'text-red-400'}>
                                        {parsedBreakdown.finalScore !== null ? parsedBreakdown.finalScore.toFixed(1) : score.toFixed(1)}
                                    </span>
                                    {parsedBreakdown.totalMax > 0 && (
                                        <span className="text-gray-600"> / {parsedBreakdown.totalMax} pts</span>
                                    )}
                                </span>
                            </div>

                            {/* Scored rows */}
                            <div className="space-y-3">
                                {parsedBreakdown.scoredRows.map(row => {
                                    const isNegative = row.contribution < 0
                                    const isDrugLikeness = row.key === 'drug_likeness'

                                    return (
                                        <div key={row.key}>
                                            {/* Row label + pts */}
                                            <div className="flex justify-between text-xs mb-1 items-center">
                                                <span className="text-gray-500 capitalize flex items-center">
                                                    {row.key.replace(/_/g, ' ')}
                                                    {BREAKDOWN_TIPS[row.key] && <Tip text={BREAKDOWN_TIPS[row.key]} />}
                                                </span>
                                                <span className={`font-mono ${isNegative
                                                        ? 'text-red-400'
                                                        : row.maxPossible > 0
                                                            ? row.pct >= 70 ? 'text-emerald-400' : row.pct >= 40 ? 'text-yellow-400' : 'text-red-400'
                                                            : 'text-gray-400'
                                                    }`}>
                                                    {row.contribution % 1 === 0
                                                        ? row.contribution.toFixed(0)
                                                        : row.contribution.toFixed(1)}
                                                    {row.maxPossible > 0
                                                        ? <span className="text-gray-600"> / {row.maxPossible} pts</span>
                                                        : <span className="text-gray-700"> pts</span>}
                                                </span>
                                            </div>

                                            {/* Progress bar */}
                                            {row.maxPossible > 0 && !isNegative && (
                                                <div className="h-1.5 bg-gray-800 rounded-full overflow-hidden mb-1.5">
                                                    <div className={`h-full rounded-full transition-all ${row.pct >= 70 ? 'bg-emerald-600' : row.pct >= 40 ? 'bg-yellow-600' : 'bg-red-700'}`} style={{ width: `${row.pct}%` }} />
                                                </div>
                                            )}
                                            {isNegative && (
                                                <div className="h-1.5 bg-gray-800 rounded-full overflow-hidden mb-1.5">
                                                    <div className="h-full rounded-full bg-red-800 transition-all" style={{ width: `${Math.min(Math.abs(row.contribution) * 2, 100)}%` }} />
                                                </div>
                                            )}

                                            {/* Raw value */}
                                            {row.raw && (
                                                isDrugLikeness
                                                    ? <DrugLikenessChips raw={row.raw} />
                                                    : <p className="text-[10px] text-gray-600 leading-relaxed">{row.raw}</p>
                                            )}
                                        </div>
                                    )
                                })}
                            </div>

                            {/* Binding prefilter gate info pill */}
                            {parsedBreakdown.gateRow && (
                                <div className="mt-3 flex items-center gap-2 px-2.5 py-2 rounded-lg border border-gray-700/60 bg-gray-800/30">
                                    <Filter className="w-3 h-3 text-gray-600 flex-shrink-0" />
                                    <div className="min-w-0">
                                        <span className="text-[10px] text-gray-600 uppercase tracking-wider font-medium">ML Binding Pre-filter</span>
                                        <span className={`ml-2 text-[10px] font-semibold ${parsedBreakdown.gateRow.passed ? 'text-emerald-500' : 'text-red-500'}`}>
                                            {parsedBreakdown.gateRow.passed ? '✓ passed' : '✗ failed'}
                                        </span>
                                        {parsedBreakdown.gateRow.raw && <p className="text-[10px] text-gray-700 mt-0.5 font-mono">{parsedBreakdown.gateRow.raw}</p>}
                                    </div>
                                    <span className="ml-auto text-[10px] text-gray-700 flex-shrink-0">gate only · 0 pts</span>
                                </div>
                            )}

                            {/* Disabled components notice */}
                            {parsedBreakdown.disabledComponents.length > 0 && (
                                <p className="text-[10px] text-gray-700 mt-2">
                                    Disabled: {parsedBreakdown.disabledComponents.join(', ')} — weight redistributed to active components.
                                </p>
                            )}

                            {/* MW fragment penalty notice */}
                            {parsedBreakdown.penaltyApplied && (
                                <p className="text-xs text-yellow-500 mt-3">⚠ Fragment penalty applied — MW &lt; 200 Da, final score halved.</p>
                            )}
                        </div>
                    )}
                </div>
            )}
        </div>
    )
}
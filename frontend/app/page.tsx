import { FlaskConical } from "lucide-react";
import InputForm from "@/components/InputForm";


const PIPELINE_STEPS = [
  { label: "Analogue Generation", desc: "RDKit BRICS" },
  { label: "Drug-likeness Filter", desc: "Lipinski RO5 + ESOL" },
  { label: "ADMET Toxicity", desc: "ADMET-AI" },
  { label: "Binding Pre-filter", desc: "DeepChem AttentiveFP" },
  { label: "Molecular Docking", desc: "GNINA 1.3" },
  { label: "Retrosynthesis", desc: "AiZynthFinder" },
];


export default function HomePage() {
  return (
    <div className="grid grid-cols-1 lg:grid-cols-5 gap-12 items-start max-w-5xl mx-auto">

      {/* ── Left: Minimal Context ───────────────────────────── */}
      <div className="lg:col-span-2 lg:sticky lg:top-24">
        <h1 className="text-3xl font-semibold tracking-tight mb-4 text-gray-100">
          Compound Screening
        </h1>

        <p className="text-gray-400 text-sm leading-relaxed mb-8">
          Submit a base SMILES and a target protein structure. The pipeline will automatically generate analogues, apply ADMET filters, perform CNN-scored docking, and assess synthetic feasibility.
        </p>

        <div className="space-y-1">
          <p className="text-[10px] font-medium text-gray-500 uppercase tracking-widest mb-3">
            Execution Flow
          </p>
          {PIPELINE_STEPS.map(({ label, desc }, i) => (
            <div key={label} className="flex items-start gap-3 py-1.5 text-sm">
              <span className="text-gray-600 font-mono text-xs mt-0.5">{i + 1}</span>
              <div>
                <p className="text-gray-300 font-medium">{label}</p>
                <p className="text-gray-500 text-xs">{desc}</p>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* ── Right: Input Form ───────────────────────────── */}
      <div className="lg:col-span-3">
        <div className="flex items-center gap-2 mb-4">
          <FlaskConical className="w-4 h-4 text-emerald-500" />
          <h2 className="text-base font-medium text-gray-200">Pipeline Parameters</h2>
        </div>
        <InputForm />
      </div>

    </div>
  );
}
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
    <div className="grid grid-cols-1 lg:grid-cols-5 gap-10 lg:gap-14 items-start max-w-5xl mx-auto">
      {/* ── Left: context ─────────────────────────────────────── */}
      <div className="lg:col-span-2 lg:sticky lg:top-24">
        <p className="eyebrow mb-3">Virtual screening</p>

        <h1 className="text-[2rem] sm:text-4xl font-semibold tracking-[-0.03em] leading-[1.1] text-gray-50 mb-4">
          Compound
          <br />
          <span className="gradient-text">Screening</span>
        </h1>

        <p className="text-gray-500 text-sm leading-relaxed mb-8 max-w-sm">
          Submit a base SMILES string and a target protein structure. The
          pipeline generates analogues, applies ADMET filters, performs
          CNN-scored docking, and assesses synthetic feasibility.
        </p>

        <p className="eyebrow mb-4">Execution flow</p>

        <ol className="flow-rail space-y-0.5">
          {PIPELINE_STEPS.map(({ label, desc }, i) => (
            <li
              key={label}
              className="flow-step flex items-start gap-3 py-2 group"
            >
              <span className="flow-marker mt-px">{i + 1}</span>
              <div className="min-w-0">
                <p className="text-sm text-gray-300 font-medium leading-tight group-hover:text-gray-100 transition-colors">
                  {label}
                </p>
                <p className="text-gray-600 text-xs mt-0.5 font-mono">{desc}</p>
              </div>
            </li>
          ))}
        </ol>
      </div>

      {/* ── Right: input form ─────────────────────────────────── */}
      <div className="lg:col-span-3">
        <div className="flex items-center gap-2 mb-4">
          <FlaskConical className="w-4 h-4 text-emerald-500" />
          <h2 className="text-sm font-semibold text-gray-200 tracking-tight">
            Pipeline Parameters
          </h2>
        </div>
        <InputForm />
      </div>
    </div>
  );
}

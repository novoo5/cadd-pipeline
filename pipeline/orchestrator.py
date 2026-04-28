import asyncio
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from backend.models.schemas import (
    JobRequest,
    JobStatusResponse,
    JobResultsResponse,
    PipelineStepInfo,
    StepStatus,
    CompoundResult,
    BindingSiteInfo,
)
from backend.pipeline.analogue_gen import generate_analogues
from backend.pipeline.drug_likeness import (
    filter_by_drug_likeness,
    filter_by_solubility,
)
from backend.pipeline.admet import run_admet_batch
from backend.pipeline.binding_prefilter import run_binding_prefilter
from backend.pipeline.protein_prep import prepare_receptor
from backend.pipeline.docking import dock_compounds, cleanup_docking_files
from backend.pipeline.retrosynthesis import run_retrosynthesis_batch
from backend.pipeline.scorer import score_and_rank
from backend.utils.pdb_fetcher import fetch_pdb_by_id, save_pdb_from_content
from backend.utils.smiles_validator import validate_and_canonicalize


JOBS_DIR = Path("storage/jobs")


class PipelineOrchestrator:
    def __init__(self, job_id: str, request: JobRequest):
        self.job_id = job_id
        self.request = request
        self.job_dir = JOBS_DIR / job_id
        self.job_dir.mkdir(parents=True, exist_ok=True)
        self.created_at = datetime.now(timezone.utc)
        self.steps: list[PipelineStepInfo] = []
        self._init_steps()

    # ── Step initialisation ───────────────────────────────────────────────────

    def _init_steps(self):
        step_names = [
            "Analogue Generation",
            "Drug-likeness Filter (Lipinski + ESOL)",
            "ADMET Toxicity Filter",
            "ML Binding Pre-filter",
            "Protein Preparation",
            "Molecular Docking (Vina)",
            "Retrosynthesis (SA Score)",
            "Scoring & Ranking",
        ]
        for name in step_names:
            self.steps.append(PipelineStepInfo(name=name, status=StepStatus.waiting))
        self._write_status("queued")

    # ── Status helpers ────────────────────────────────────────────────────────

    def _write_status(self, status: str, error: Optional[str] = None):
        payload = JobStatusResponse(
            job_id=self.job_id,
            status=status,
            created_at=self.created_at,
            updated_at=datetime.now(timezone.utc),
            steps=self.steps,
            error=error,
        )
        (self.job_dir / "status.json").write_text(payload.model_dump_json(indent=2))

    def _step_start(self, index: int):
        self.steps[index].status = StepStatus.running
        self._write_status("running")

    def _step_done(self, index: int, message: str, elapsed: float):
        self.steps[index].status = StepStatus.done
        self.steps[index].message = message
        self.steps[index].duration_seconds = round(elapsed, 2)
        self.steps[index].progress_current = None
        self.steps[index].progress_total = None
        self._write_status("running")

    def _step_skip(self, index: int, reason: str = "Disabled by user"):
        self.steps[index].status = StepStatus.skipped
        self.steps[index].message = reason
        self._write_status("running")

    def _step_fail(self, index: int, error: str):
        self.steps[index].status = StepStatus.failed
        self.steps[index].message = error
        self._write_status("failed", error=error)

    def _update_step_progress(self, index: int, current: int, total: int):
        self.steps[index].progress_current = current
        self.steps[index].progress_total = total
        self._write_status("running")

    # ── Main pipeline ─────────────────────────────────────────────────────────

    async def run(self):
        req = self.request
        compound_results: dict[str, CompoundResult] = {}
        binding_site_info: Optional[BindingSiteInfo] = None

        locked_scaffold = getattr(req, "locked_scaffold_smarts", None) or None

        # ── Routing flags ─────────────────────────────────────────────────────
        tox_mode = req.toxicity_report_only
        direct_only = req.direct_score_only and not req.toxicity_report_only

        # ── admet_config — use what the frontend sent, or None for defaults ──
        admet_config = getattr(req, "admet_config", None)

        try:
            is_valid, canonical_smiles, err = validate_and_canonicalize(req.smiles)
            if not is_valid:
                raise ValueError(f"Invalid SMILES: {err}")

            # ── STEP 1: Analogue Generation ───────────────────────────────────
            t = time.time()
            self._step_start(0)

            def analogue_progress(current: int, total: int):
                self._update_step_progress(0, current, total)

            if tox_mode:
                if req.direct_score_only:
                    analogues = [canonical_smiles]
                    self._step_done(
                        0,
                        "Direct Score + Toxicity Report — base compound only queued for ADMET (no analogues)",
                        time.time() - t,
                    )
                elif req.num_analogues > 0:
                    analogues = await asyncio.to_thread(
                        generate_analogues,
                        canonical_smiles,
                        req.num_analogues,
                        req.solubility_filter.value,
                        on_progress=analogue_progress,
                        locked_scaffold_smarts=locked_scaffold,
                    )
                    if not analogues:
                        analogues = [canonical_smiles]
                    if canonical_smiles not in analogues:
                        analogues = [canonical_smiles] + analogues

                    scaffold_note = " · scaffold locked" if locked_scaffold else ""
                    self._step_done(
                        0,
                        (
                            f"Toxicity Report mode — {len(analogues)} compound(s) "
                            f"queued for ADMET screening "
                            f"(base + {len(analogues) - 1} analogues){scaffold_note}"
                        ),
                        time.time() - t,
                    )
                else:
                    analogues = [canonical_smiles]
                    self._step_done(
                        0,
                        "Toxicity Report mode — base compound only (no analogues requested)",
                        time.time() - t,
                    )

            elif direct_only:
                analogues = [canonical_smiles]
                self._step_done(
                    0,
                    "Direct scoring mode — input compound scored as-is (no analogues generated)",
                    time.time() - t,
                )

            else:
                analogues = await asyncio.to_thread(
                    generate_analogues,
                    canonical_smiles,
                    req.num_analogues,
                    req.solubility_filter.value,
                    on_progress=analogue_progress,
                    locked_scaffold_smarts=locked_scaffold,
                )
                if not analogues:
                    raise RuntimeError("Analogue generation produced 0 compounds.")

                sol_note = (
                    f" — solubility filter: {req.solubility_filter.value}"
                    if req.solubility_filter.value != "all"
                    else ""
                )
                scaffold_note = " · scaffold locked" if locked_scaffold else ""
                self._step_done(
                    0,
                    f"{len(analogues)}/{req.num_analogues} analogues generated{sol_note}{scaffold_note}",
                    time.time() - t,
                )

            for smi in analogues:
                compound_results[smi] = CompoundResult(smiles=smi, canonical_smiles=smi)

            # ── STEP 2: Drug-likeness Filter ──────────────────────────────────
            t = time.time()
            self._step_start(1)

            if tox_mode:
                self._step_skip(1, "Skipped — toxicity report mode screens all compounds")
                surviving = analogues

            elif req.pipeline_steps.drug_likeness:
                passed_lipinski, lipinski_results = filter_by_drug_likeness(
                    analogues,
                    mw_min=req.mw_min,
                    mw_max=req.mw_max,
                    max_violations=req.max_lipinski_violations,
                )
                for smi, result in zip(analogues, lipinski_results):
                    compound_results[smi].lipinski = result

                sol_filter = req.solubility_filter.value if req.solubility_filter else "all"
                if sol_filter != "all":
                    before_sol = len(passed_lipinski)
                    passed_lipinski = filter_by_solubility(
                        passed_lipinski, compound_results, sol_filter
                    )
                    sol_dropped = before_sol - len(passed_lipinski)
                    sol_label = f" · {sol_dropped} removed by solubility filter ({sol_filter})"
                else:
                    sol_label = ""

                if req.max_lipinski_violations is None:
                    filter_label = "Lipinski filter ignored (all compounds pass)"
                else:
                    filter_label = f"≤{req.max_lipinski_violations} violation(s) allowed"

                self._step_done(
                    1,
                    f"{len(passed_lipinski)}/{len(analogues)} passed Lipinski + ESOL ({filter_label}){sol_label}",
                    time.time() - t,
                )
                surviving = passed_lipinski
            else:
                self._step_skip(1)
                surviving = analogues

            if not surviving and not tox_mode:
                raise RuntimeError("No compounds survived the drug-likeness filter.")

            # ── STEP 3: ADMET ─────────────────────────────────────────────────
            t_pdb = time.time()

            if tox_mode:
                self._step_start(2)
                t_admet = time.time()

                # ← admet_config forwarded from request
                passed_admet, admet_results = await asyncio.to_thread(
                    run_admet_batch, surviving, admet_config
                )
                for smi, result in zip(surviving, admet_results):
                    compound_results[smi].admet = result

                total_flags = sum(len(r.flags) for r in admet_results)
                high_flags = sum(
                    1 for r in admet_results for f in r.flags if f.severity == "high"
                )
                self._step_done(
                    2,
                    (
                        f"{len(passed_admet)}/{len(surviving)} passed ADMET — "
                        f"{total_flags} flag(s) raised ({high_flags} high-severity)"
                    ),
                    time.time() - t_admet,
                )

                self._step_skip(3, "Skipped — toxicity report mode")
                self._step_skip(4, "Skipped — no docking in toxicity report mode")
                self._step_skip(5, "Skipped — toxicity report mode")
                self._step_skip(6, "Skipped — toxicity report mode")

                surviving = passed_admet
                receptor_pdbqt = None

            else:
                self._step_start(4)
                if req.pdb_id:
                    raw_pdb_path = await fetch_pdb_by_id(req.pdb_id, self.job_id)
                else:
                    raw_pdb_path = save_pdb_from_content(req.pdb_content, self.job_id)

                if req.pipeline_steps.admet:
                    self._step_start(2)
                    t_admet = time.time()

                    # ← admet_config forwarded from request
                    admet_task = asyncio.to_thread(
                        run_admet_batch, surviving, admet_config
                    )
                    protprep_task = asyncio.to_thread(
                        prepare_receptor, raw_pdb_path, self.job_id
                    )
                    (passed_admet, admet_results), receptor_pdbqt = await asyncio.gather(
                        admet_task, protprep_task
                    )
                    for smi, result in zip(surviving, admet_results):
                        compound_results[smi].admet = result

                    total_flags = sum(len(r.flags) for r in admet_results)
                    high_flags = sum(
                        1 for r in admet_results for f in r.flags if f.severity == "high"
                    )
                    self._step_done(
                        2,
                        (
                            f"{len(passed_admet)}/{len(surviving)} passed ADMET — "
                            f"{total_flags} flag(s) raised ({high_flags} high-severity)"
                        ),
                        time.time() - t_admet,
                    )
                    self._step_done(
                        4,
                        f"Receptor prepared → {Path(receptor_pdbqt).name}",
                        time.time() - t_pdb,
                    )
                    surviving = passed_admet
                else:
                    self._step_skip(2)
                    receptor_pdbqt = await asyncio.to_thread(
                        prepare_receptor, raw_pdb_path, self.job_id
                    )
                    self._step_done(
                        4,
                        f"Receptor prepared → {Path(receptor_pdbqt).name}",
                        time.time() - t_pdb,
                    )

                if not surviving:
                    raise RuntimeError("No compounds survived the ADMET filter.")

                # ── STEP 4: ML Binding Pre-filter ─────────────────────────────
                t = time.time()
                self._step_start(3)
                if req.pipeline_steps.binding_prefilter:
                    top_smiles, prefilter_results = await asyncio.to_thread(
                        run_binding_prefilter, surviving, req.max_docking_compounds
                    )
                    for smi, result in zip(surviving, prefilter_results):
                        compound_results[smi].binding_prefilter = result
                    self._step_done(
                        3,
                        f"Top {len(top_smiles)}/{len(surviving)} selected for docking",
                        time.time() - t,
                    )
                    surviving = top_smiles
                else:
                    self._step_skip(3)
                    surviving = surviving[: req.max_docking_compounds]

                if not surviving:
                    raise RuntimeError("No compounds survived the binding pre-filter.")

                # ── STEP 5: Vina Docking ──────────────────────────────────────
                t = time.time()
                self._step_start(5)
                if req.pipeline_steps.docking:
                    def docking_progress(current: int, total: int):
                        self._update_step_progress(5, current, total)

                    passed_docking, docking_results, binding_site_info = await asyncio.to_thread(
                        dock_compounds,
                        surviving,
                        receptor_pdbqt,
                        self.job_id,
                        raw_pdb_path,
                        req.binding_site_mode,
                        req.docking_speed,
                        req.binding_site_coords,
                        req.binding_site_residues,
                        on_progress=docking_progress,
                    )
                    for smi, result in zip(surviving, docking_results):
                        compound_results[smi].docking = result
                    self._step_done(
                        5,
                        f"{len(passed_docking)}/{len(surviving)} passed docking (affinity ≤ -5.0 kcal/mol)",
                        time.time() - t,
                    )
                    surviving = passed_docking if passed_docking else surviving[:3]
                else:
                    self._step_skip(5)

                # ── STEP 6: Retrosynthesis ────────────────────────────────────
                t = time.time()
                self._step_start(6)
                if req.pipeline_steps.retrosynthesis:
                    retro_results = await run_retrosynthesis_batch(surviving)
                    for smi, result in zip(surviving, retro_results):
                        compound_results[smi].retrosynthesis = result
                    feasible_count = sum(1 for r in retro_results if r.feasible)
                    self._step_done(
                        6,
                        f"{feasible_count}/{len(surviving)} compounds have feasible synthesis routes",
                        time.time() - t,
                    )
                else:
                    self._step_skip(6)

            # ── STEP 7: Scoring & Ranking (always runs) ───────────────────────
            t = time.time()
            self._step_start(7)
            ranked = score_and_rank(
                list(compound_results.values()),
                pipeline_steps=req.pipeline_steps,
            )
            self._step_done(7, f"Ranked {len(ranked)} compounds", time.time() - t)

            self._write_results(ranked, canonical_smiles, len(analogues), binding_site_info)
            self._write_status("done")

            if not tox_mode and req.pipeline_steps.docking:
                cleanup_docking_files(self.job_id)

        except Exception as e:
            print(f"[Orchestrator] Job {self.job_id} failed: {e}")
            failed_step = next(
                (i for i, s in enumerate(self.steps) if s.status == StepStatus.running),
                7,
            )
            self._step_fail(failed_step, str(e))
            self._write_status("failed", error=str(e))
            raise

    # ── Results writer ────────────────────────────────────────────────────────

    def _write_results(
        self,
        ranked: list[CompoundResult],
        base_smiles: str,
        total_generated: int,
        binding_site_info: Optional[BindingSiteInfo] = None,
    ):
        def _count(attr: str) -> int:
            return sum(
                1
                for c in ranked
                if getattr(c, attr) and getattr(getattr(c, attr), "passed", False)
            )

        results = JobResultsResponse(
            job_id=self.job_id,
            base_smiles=base_smiles,
            total_analogues_generated=total_generated,
            compounds_after_lipinski=_count("lipinski"),
            compounds_after_admet=_count("admet"),
            compounds_after_prefilter=_count("binding_prefilter"),
            compounds_docked=_count("docking"),
            final_ranked_compounds=ranked,
            binding_site_info=binding_site_info,
            created_at=self.created_at,
            completed_at=datetime.now(timezone.utc),
        )
        (self.job_dir / "results.json").write_text(results.model_dump_json(indent=2))
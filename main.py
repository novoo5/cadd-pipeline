import uuid
import json
from pathlib import Path
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from typing import Optional

from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from backend.models.schemas import JobRequest, JobSubmitResponse
from backend.routers import jobs, health
from backend.worker import enqueue_job


JOBS_DIR = Path("storage/jobs")
JOBS_DIR.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://novoo5-cadd-backend.hf.space"


# ── Initial status writer ─────────────────────────────────────────────────────

def _write_initial_status(job_id: str, request: JobRequest):
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    now = datetime.now(timezone.utc).isoformat()
    ps  = request.pipeline_steps

    # In toxicity_report_only mode steps 2/4/5/6 are pre-marked skipped
    tox_only = request.toxicity_report_only

    steps = [
        {"name": "Analogue Generation",               "status": "waiting"},
        {"name": "Drug-likeness Filter (Lipinski + ESOL)",
         "status": "skipped" if tox_only else ("waiting" if ps.drug_likeness else "skipped")},
        {"name": "ADMET Toxicity Filter",              "status": "waiting" if ps.admet else "skipped"},
        {"name": "ML Binding Pre-filter",
         "status": "skipped" if tox_only else ("waiting" if ps.binding_prefilter else "skipped")},
        {"name": "Protein Preparation",
         "status": "skipped" if tox_only else ("waiting" if ps.docking else "skipped")},
        {"name": "Molecular Docking (Vina)",
         "status": "skipped" if tox_only else ("waiting" if ps.docking else "skipped")},
        {"name": "Retrosynthesis (SA Score)",
         "status": "skipped" if tox_only else ("waiting" if ps.retrosynthesis else "skipped")},
        {"name": "Scoring & Ranking",                  "status": "waiting"},
    ]

    status_payload = {
        "job_id":     job_id,
        "status":     "queued",
        "created_at": now,
        "updated_at": now,
        "steps":      steps,
        "error":      None,
    }

    (job_dir / "status.json").write_text(json.dumps(status_payload, indent=2))
    (job_dir / "request.json").write_text(request.model_dump_json(indent=2))


# ── Lifespan ──────────────────────────────────────────────────────────────────

@asynccontextmanager
async def lifespan(app: FastAPI):
    print("[CADD] Backend starting up...")
    yield
    print("[CADD] Backend shutting down...")


# ── App ───────────────────────────────────────────────────────────────────────

app = FastAPI(
    title="CADD Pipeline API",
    description=(
        "Computer-Aided Drug Design pipeline: analogue generation → "
        "drug-likeness → ADMET → binding pre-filter → Vina docking → "
        "retrosynthesis → scoring."
    ),
    version="1.1.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://*.vercel.app",
        "http://localhost:3000",
        "http://192.168.1.4:3000",
        "http://127.0.0.1:3000",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(jobs.router,   prefix="/api/v1", tags=["Jobs"])
app.include_router(health.router, prefix="/api/v1", tags=["Health"])


# ── Root ──────────────────────────────────────────────────────────────────────

@app.get("/")
async def root():
    return {
        "name":    "CADD Pipeline API",
        "version": "1.1.0",
        "status":  "running",
        "docs":    "/docs",
    }


# ── JSON submission ───────────────────────────────────────────────────────────

@app.post("/api/v1/jobs", response_model=JobSubmitResponse)
async def submit_job(request: JobRequest):
    """
    Submit a pipeline job via JSON body.

    New fields:
      solubility_filter    — "soluble_only" | "allow_slightly" | "all"
      toxicity_report_only — true = ADMET-only, no PDB required
    """
    job_id = str(uuid.uuid4())
    _write_initial_status(job_id, request)
    await enqueue_job(job_id, request.model_dump())

    return JobSubmitResponse(
        job_id      = job_id,
        message     = "Job queued successfully. Use the status URL to track progress.",
        status_url  = f"{BASE_URL}/api/v1/jobs/{job_id}/status",
        results_url = f"{BASE_URL}/api/v1/jobs/{job_id}/results",
    )


# ── File-upload submission ────────────────────────────────────────────────────

@app.post("/api/v1/jobs/upload", response_model=JobSubmitResponse)
async def submit_job_with_file(
    smiles:                  str            = Form(...),
    num_analogues:           int            = Form(25),       # 0–1000
    docking_speed:           str            = Form("balanced"),
    binding_site_mode:       str            = Form("auto"),
    direct_score_only:       bool           = Form(False),
    mw_min:                  float          = Form(200.0),
    mw_max:                  float          = Form(500.0),
    # Frontend sends -1 to represent None ("ignore filter")
    max_lipinski_violations: Optional[int]  = Form(1),
    # ── NEW ──────────────────────────────────────────────────────────────────
    solubility_filter:       str            = Form("all"),
    toxicity_report_only:    bool           = Form(False),
    # ─────────────────────────────────────────────────────────────────────────
    pdb_id:                  Optional[str]  = Form(None),
    pdb_file:                UploadFile     = File(None),
):
    """
    Submit a pipeline job via multipart/form-data (for .pdb file uploads).

    solubility_filter:
        "soluble_only"    — only logS > -3 analogues are generated
        "allow_slightly"  — logS > -5 (also slightly soluble)
        "all"             — no filter (default)

    toxicity_report_only:
        true  — skip docking/retrosynthesis; run ADMET only; PDB not required
        false — full pipeline (default)
    """
    from backend.models.schemas import (
        PipelineSteps, BindingSiteMode, DockingSpeed, SolubilityFilterMode
    )

    # PDB required only when docking is going to run
    needs_pdb = not toxicity_report_only
    if needs_pdb and pdb_id is None and pdb_file is None:
        raise HTTPException(
            status_code=422,
            detail=(
                "Either pdb_id or pdb_file must be provided. "
                "Not required when toxicity_report_only=true."
            ),
        )

    # Parse uploaded PDB file
    pdb_content: Optional[str] = None
    if pdb_file is not None:
        if not pdb_file.filename.endswith(".pdb"):
            raise HTTPException(
                status_code=422,
                detail="Uploaded file must have a .pdb extension.",
            )
        raw_bytes   = await pdb_file.read()
        pdb_content = raw_bytes.decode("utf-8", errors="ignore")
        if len(pdb_content) < 100:
            raise HTTPException(
                status_code=422,
                detail="Uploaded PDB file appears to be empty or too small.",
            )

    # Validate solubility_filter value before passing to schema
    try:
        sol_filter = SolubilityFilterMode(solubility_filter)
    except ValueError:
        raise HTTPException(
            status_code=422,
            detail=f"Invalid solubility_filter '{solubility_filter}'. "
                   f"Must be one of: soluble_only, allow_slightly, all.",
        )

    # -1 from frontend = None ("ignore Lipinski filter entirely")
    resolved_violations = None if max_lipinski_violations == -1 else max_lipinski_violations

    request = JobRequest(
        smiles                  = smiles,
        pdb_id                  = pdb_id,
        pdb_content             = pdb_content,
        num_analogues           = 0 if direct_score_only else num_analogues,
        direct_score_only       = direct_score_only,
        solubility_filter       = sol_filter,           # ← NEW
        toxicity_report_only    = toxicity_report_only, # ← NEW
        pipeline_steps          = PipelineSteps(),
        binding_site_mode       = BindingSiteMode(binding_site_mode),
        docking_speed           = DockingSpeed(docking_speed),
        mw_min                  = mw_min,
        mw_max                  = mw_max,
        max_lipinski_violations = resolved_violations,
    )

    job_id = str(uuid.uuid4())
    _write_initial_status(job_id, request)
    await enqueue_job(job_id, request.model_dump())

    return JobSubmitResponse(
        job_id      = job_id,
        message     = "Job queued from file upload successfully.",
        status_url  = f"{BASE_URL}/api/v1/jobs/{job_id}/status",
        results_url = f"{BASE_URL}/api/v1/jobs/{job_id}/results",
    )


# ── Global exception handler ──────────────────────────────────────────────────

@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    return JSONResponse(
        status_code=500,
        content={
            "error":  "Internal server error",
            "detail": str(exc),
        },
    )

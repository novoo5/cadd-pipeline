import uuid
import json
import time
from pathlib import Path
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from backend.models.schemas import JobRequest, JobSubmitResponse
from backend.routers import jobs, health
from backend.worker import enqueue_job


JOBS_DIR = Path("storage/jobs")
JOBS_DIR.mkdir(parents=True, exist_ok=True)


# ── Helper: write queued status.json before ARQ picks it up ──────────────────

def _write_initial_status(job_id: str, request: JobRequest):
    """
    Writes status.json immediately at submission time so /status never 404s
    while the job is sitting in the ARQ queue waiting for a worker slot.
    Without this, jobs returning 404 when queue has 3+ jobs ahead.
    """
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    now = datetime.now(timezone.utc).isoformat()

    # Build the step list — mark skipped steps based on pipeline_steps config
    ps = request.pipeline_steps
    steps = [
        {"name": "Drug-likeness Filter", "status": "waiting"},
        {"name": "ADMET Toxicity Filter", "status": "waiting" if ps.admet else "skipped"},
        {"name": "ML Binding Pre-filter", "status": "waiting" if ps.binding_prefilter else "skipped"},
        {"name": "Protein Preparation",   "status": "waiting" if ps.docking else "skipped"},
        {"name": "Molecular Docking",     "status": "waiting" if ps.docking else "skipped"},
        {"name": "Retrosynthesis",        "status": "waiting" if ps.retrosynthesis else "skipped"},
        {"name": "Scoring & Ranking",     "status": "waiting"},
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


app = FastAPI(
    title="CADD Pipeline API",
    description=(
        "Computer-Aided Drug Design pipeline: analogue generation → "
        "drug-likeness → ADMET → binding pre-filter → GNINA docking → "
        "ASKCOS retrosynthesis → scoring."
    ),
    version="1.0.0",
    lifespan=lifespan,
)


# ── CORS ──────────────────────────────────────────────────────────────────────
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


# ── Routers ───────────────────────────────────────────────────────────────────
app.include_router(jobs.router,   prefix="/api/v1", tags=["Jobs"])
app.include_router(health.router, prefix="/api/v1", tags=["Health"])


# ── Root ──────────────────────────────────────────────────────────────────────
@app.get("/")
async def root():
    return {
        "name": "CADD Pipeline API",
        "version": "1.0.0",
        "status": "running",
        "docs": "/docs",
    }


# ── Job submission via JSON ───────────────────────────────────────────────────
@app.post("/api/v1/jobs", response_model=JobSubmitResponse)
async def submit_job(request: JobRequest):
    """
    Submit a new pipeline job via JSON body.
    Pass "direct_score_only": true to skip analogue generation.
    Pass "mw_min" and "mw_max" to override the MW filter range (default 200–500 Da).
    """
    job_id = str(uuid.uuid4())

    # ✅ Write status.json BEFORE enqueuing — fixes 404 race condition
    # Previously: only request.json was written here, status.json was created
    # by the orchestrator only after the worker picked it up. With a full queue
    # (3+ jobs ahead), the frontend would hit /status, get 404, and give up.
    _write_initial_status(job_id, request)

    await enqueue_job(job_id, request.model_dump())

    base_url = "https://novoo5-cadd-backend.hf.space"
    return JobSubmitResponse(
        job_id=job_id,
        message="Job queued successfully. Use the status URL to track progress.",
        status_url=f"{base_url}/api/v1/jobs/{job_id}/status",
        results_url=f"{base_url}/api/v1/jobs/{job_id}/results",
    )


# ── Job submission via file upload (multipart form) ───────────────────────────
@app.post("/api/v1/jobs/upload", response_model=JobSubmitResponse)
async def submit_job_with_file(
    smiles: str = Form(...),
    num_analogues: int = Form(25),
    docking_speed: str = Form("balanced"),
    binding_site_mode: str = Form("auto"),
    direct_score_only: bool = Form(False),
    mw_min: float = Form(200.0),
    mw_max: float = Form(500.0),
    pdb_id: str = Form(None),
    pdb_file: UploadFile = File(None),
):
    """
    Alternative submission for .pdb file uploads (multipart/form-data).
    Set direct_score_only=true to bypass analogue generation.
    Set mw_min/mw_max to override MW filter range.
    """
    if pdb_id is None and pdb_file is None:
        raise HTTPException(
            status_code=422,
            detail="Either pdb_id or pdb_file must be provided.",
        )

    pdb_content = None
    if pdb_file is not None:
        if not pdb_file.filename.endswith(".pdb"):
            raise HTTPException(
                status_code=422,
                detail="Uploaded file must have a .pdb extension.",
            )
        raw_bytes = await pdb_file.read()
        pdb_content = raw_bytes.decode("utf-8", errors="ignore")
        if len(pdb_content) < 100:
            raise HTTPException(
                status_code=422,
                detail="Uploaded PDB file appears to be empty or too small.",
            )

    from backend.models.schemas import PipelineSteps, BindingSiteMode, DockingSpeed

    request = JobRequest(
        smiles=smiles,
        pdb_id=pdb_id,
        pdb_content=pdb_content,
        num_analogues=0 if direct_score_only else num_analogues,
        direct_score_only=direct_score_only,
        pipeline_steps=PipelineSteps(),
        binding_site_mode=BindingSiteMode(binding_site_mode),
        docking_speed=DockingSpeed(docking_speed),
        mw_min=mw_min,
        mw_max=mw_max,
    )

    job_id = str(uuid.uuid4())

    # ✅ Write status.json BEFORE enqueuing — same fix as JSON endpoint
    _write_initial_status(job_id, request)

    await enqueue_job(job_id, request.model_dump())

    base_url = "https://novoo5-cadd-backend.hf.space"
    return JobSubmitResponse(
        job_id=job_id,
        message="Job queued from file upload successfully.",
        status_url=f"{base_url}/api/v1/jobs/{job_id}/status",
        results_url=f"{base_url}/api/v1/jobs/{job_id}/results",
    )


# ── Global exception handler ──────────────────────────────────────────────────
@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    return JSONResponse(
        status_code=500,
        content={
            "error": "Internal server error",
            "detail": str(exc),
        },
    )

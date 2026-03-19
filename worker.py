import json
import asyncio
from pathlib import Path
from typing import Any

# ARQ — async Redis Queue for Python
# pip install arq
# Redis must be running (added to Dockerfile)
from arq import create_pool
from arq.connections import RedisSettings

from backend.models.schemas import JobRequest
from backend.pipeline.orchestrator import PipelineOrchestrator

JOBS_DIR = Path("storage/jobs")

# Redis connection — runs on localhost inside the Docker container
REDIS_SETTINGS = RedisSettings(host="localhost", port=6379)


# ── Job function (executed by the ARQ worker process) ────────────────────────

async def run_pipeline_job(ctx: dict, job_id: str, request_data: dict):
    """
    ARQ worker function — this is what actually runs in the background.
    ARQ calls this with the serialized job_id + request_data dict.
    """
    print(f"[Worker] Starting job {job_id}")
    try:
        request = JobRequest(**request_data)
        orchestrator = PipelineOrchestrator(job_id=job_id, request=request)
        await orchestrator.run()
        print(f"[Worker] Job {job_id} completed successfully.")
    except Exception as e:
        print(f"[Worker] Job {job_id} failed: {e}")
        # Status is already written to disk by the orchestrator
        # No re-raise needed — ARQ will mark the job as failed automatically


# ── Enqueue helper (called from FastAPI endpoints) ───────────────────────────

async def enqueue_job(job_id: str, request_data: dict):
    """
    Pushes a new pipeline job into the Redis queue.
    FastAPI calls this immediately after receiving a submission request.
    Returns instantly — the worker picks it up asynchronously.
    """
    pool = await create_pool(REDIS_SETTINGS)
    await pool.enqueue_job(
        "run_pipeline_job",
        job_id,
        request_data,
        _job_id=job_id,         # use our own job_id as the ARQ job key
        _queue_name="cadd",
    )
    await pool.aclose()


# ── ARQ WorkerSettings (used when running: python -m arq backend.worker.WorkerSettings) ──

class WorkerSettings:
    """
    ARQ reads this class to configure the worker process.
    Start the worker with:
        python -m arq backend.worker.WorkerSettings
    """
    functions = [run_pipeline_job]
    redis_settings = REDIS_SETTINGS
    queue_name = "cadd"
    max_jobs = 1            # one job at a time — free HF CPU cannot handle parallelism
    job_timeout = 3600      # 1 hour max per job (accounts for thorough docking mode)
    keep_result = 86400     # keep job result in Redis for 24 hours
    retry_jobs = False      # don't retry failed jobs — user should resubmit
    health_check_interval = 30

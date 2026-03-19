import json
from pathlib import Path
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse

from backend.models.schemas import JobStatusResponse, JobResultsResponse

router = APIRouter()

JOBS_DIR = Path("storage/jobs")


@router.get("/jobs/{job_id}/status", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    """
    Returns the current status and step-by-step progress of a pipeline job.
    The frontend polls this endpoint every 8 seconds.
    Status is read directly from disk (status.json written by the orchestrator).
    """
    _validate_job_id(job_id)

    status_path = JOBS_DIR / job_id / "status.json"
    if not status_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Job '{job_id}' not found. It may have expired or never existed.",
        )

    try:
        data = json.loads(status_path.read_text())
        return JobStatusResponse(**data)
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to read status for job '{job_id}': {e}",
        )


@router.get("/jobs/{job_id}/results", response_model=JobResultsResponse)
async def get_job_results(job_id: str):
    """
    Returns the full ranked results for a completed pipeline job.
    Only available once the job status is 'done'.
    Returns 202 Accepted if the job is still running.
    """
    _validate_job_id(job_id)

    status_path  = JOBS_DIR / job_id / "status.json"
    results_path = JOBS_DIR / job_id / "results.json"

    if not status_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Job '{job_id}' not found.",
        )

    status_data = json.loads(status_path.read_text())
    job_status  = status_data.get("status")

    if job_status == "failed":
        raise HTTPException(
            status_code=500,
            detail=(
                f"Job '{job_id}' failed during pipeline execution. "
                f"Error: {status_data.get('error', 'Unknown error')}. "
                "Please resubmit with different parameters."
            ),
        )

    if job_status in ("queued", "running"):
        # 202 = job accepted but not yet complete
        return JSONResponse(
            status_code=202,
            content={
                "message": f"Job '{job_id}' is still {job_status}. Poll /status for progress.",
                "status": job_status,
            },
        )

    if not results_path.exists():
        raise HTTPException(
            status_code=404,
            detail=(
                f"Job '{job_id}' status is 'done' but results file is missing. "
                "This is an internal error — please report it."
            ),
        )

    try:
        data = json.loads(results_path.read_text())
        return JobResultsResponse(**data)
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to parse results for job '{job_id}': {e}",
        )


@router.get("/jobs/{job_id}/score-breakdown/{compound_index}")
async def get_score_breakdown(job_id: str, compound_index: int):
    """
    Returns the detailed score breakdown for a specific compound
    in the results (by its rank index, 0-based).
    Used by the frontend expandable row.
    """
    _validate_job_id(job_id)

    results_path = JOBS_DIR / job_id / "results.json"
    if not results_path.exists():
        raise HTTPException(status_code=404, detail="Results not available yet.")

    from backend.models.schemas import CompoundResult
    from backend.pipeline.scorer import get_score_breakdown

    data = json.loads(results_path.read_text())
    compounds = data.get("final_ranked_compounds", [])

    if compound_index >= len(compounds) or compound_index < 0:
        raise HTTPException(
            status_code=404,
            detail=f"Compound index {compound_index} out of range (0–{len(compounds)-1}).",
        )

    compound = CompoundResult(**compounds[compound_index])
    breakdown = get_score_breakdown(compound)
    return JSONResponse(content=breakdown)


@router.delete("/jobs/{job_id}")
async def delete_job(job_id: str):
    """
    Deletes all files associated with a job (status, results, receptor, ligands).
    Use this to free up HF Space storage manually.
    """
    import shutil
    _validate_job_id(job_id)

    job_dir = JOBS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' not found.")

    shutil.rmtree(job_dir)
    return {"message": f"Job '{job_id}' and all associated files deleted."}


@router.get("/jobs")
async def list_recent_jobs():
    """
    Lists all jobs currently on disk (for debugging / admin purposes).
    Returns job_id, status, and created_at for each.
    """
    jobs = []
    if not JOBS_DIR.exists():
        return {"jobs": []}

    for job_dir in sorted(JOBS_DIR.iterdir(), reverse=True):
        if not job_dir.is_dir():
            continue
        status_path = job_dir / "status.json"
        if not status_path.exists():
            continue
        try:
            data = json.loads(status_path.read_text())
            jobs.append({
                "job_id":     data.get("job_id"),
                "status":     data.get("status"),
                "created_at": data.get("created_at"),
                "updated_at": data.get("updated_at"),
            })
        except Exception:
            continue

    return {"jobs": jobs[:50]}   # cap at 50


def _validate_job_id(job_id: str):
    """
    Basic UUID format validation to prevent path traversal attacks.
    A valid job_id is a UUID4 string (36 chars, 4 hyphens).
    """
    import re
    uuid4_pattern = re.compile(
        r"^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$",
        re.IGNORECASE,
    )
    if not uuid4_pattern.match(job_id):
        raise HTTPException(
            status_code=400,
            detail="Invalid job_id format. Must be a valid UUID4 string.",
        )

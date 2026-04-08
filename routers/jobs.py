import json
import re
from pathlib import Path
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse, FileResponse


from backend.models.schemas import JobStatusResponse, JobResultsResponse


router = APIRouter()


JOBS_DIR = Path("storage/jobs")


# ── File type → filename mapping (must stay in sync with docking.py) ──────────
_DOCKING_FILE_MAP = {
    "pdbqt":       ("out_{idx}.pdbqt",                "application/octet-stream"),
    "pose_pdb":    ("best_pose_{idx}_H.pdb",           "chemical/x-pdb"),
    "complex_pdb": ("complex_{idx}_PLIP_ready.pdb",    "chemical/x-pdb"),
}



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



# ── Docking file downloads ─────────────────────────────────────────────────────
#
# Three files are generated per docked compound by docking.py:
#   out_{idx}.pdbqt                — raw AutoDock Vina output (all poses)
#   best_pose_{idx}_H.pdb          — pose 1 only, with hydrogens added via obabel
#   complex_{idx}_PLIP_ready.pdb   — receptor ATOM lines + docked ligand HETATM,
#                                    merged and ready for PLIP / PyMOL
#
# All three live in:  storage/jobs/{job_id}/docking_outputs/


@router.get("/jobs/{job_id}/compounds/{compound_index}/download/{file_type}")
async def download_docking_file(job_id: str, compound_index: int, file_type: str):
    """
    Download a docking output file for a specific compound.

    file_type must be one of:
      pdbqt       → raw Vina output (.pdbqt, all poses)
      pose_pdb    → best pose with hydrogens (.pdb, ready for PLIP/PyMOL)
      complex_pdb → receptor + ligand complex (.pdb, ready for PLIP)

    Returns 404 if docking was not run or the compound index is invalid.
    Returns 503 if obabel was not installed (for pose_pdb / complex_pdb).
    """
    _validate_job_id(job_id)

    if file_type not in _DOCKING_FILE_MAP:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Unknown file_type '{file_type}'. "
                "Must be one of: pdbqt, pose_pdb, complex_pdb."
            ),
        )

    if compound_index < 0:
        raise HTTPException(status_code=400, detail="compound_index must be >= 0.")

    filename_template, media_type = _DOCKING_FILE_MAP[file_type]
    filename  = filename_template.replace("{idx}", str(compound_index))
    file_path = JOBS_DIR / job_id / "docking_outputs" / filename

    # ── Job must exist ────────────────────────────────────────────────────────
    if not (JOBS_DIR / job_id).exists():
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' not found.")

    # ── Docking output directory must exist ───────────────────────────────────
    out_dir = JOBS_DIR / job_id / "docking_outputs"
    if not out_dir.exists():
        raise HTTPException(
            status_code=404,
            detail=(
                f"No docking outputs found for job '{job_id}'. "
                "Docking may not have been run for this job, or files were cleaned up."
            ),
        )

    # ── Raw PDBQT must exist before anything else is checked ─────────────────
    if not (out_dir / f"out_{compound_index}.pdbqt").exists():
        raise HTTPException(
            status_code=404,
            detail=(
                f"Docking output not found for compound index {compound_index}. "
                "This compound may not have been docked (failed earlier filter, "
                "or index is out of range)."
            ),
        )

    # ── For pose_pdb / complex_pdb, obabel must have run successfully ─────────
    if file_type in ("pose_pdb", "complex_pdb") and not file_path.exists():
        raise HTTPException(
            status_code=503,
            detail=(
                f"File '{filename}' was not generated for compound {compound_index}. "
                "This usually means obabel is not installed on the backend. "
                "You can still download the raw .pdbqt and convert it manually with: "
                "obabel out.pdbqt -O pose.pdb -f 1 -l 1 -h"
            ),
        )

    if not file_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"File '{filename}' not found for compound {compound_index}.",
        )

    # ── Serve the file ────────────────────────────────────────────────────────
    download_name = f"compound_{compound_index}_{filename}"
    return FileResponse(
        path=str(file_path),
        media_type=media_type,
        filename=download_name,
        headers={
            "Content-Disposition": f'attachment; filename="{download_name}"',
            "Cache-Control": "no-cache",
        },
    )



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

    return {"jobs": jobs[:50]}



def _validate_job_id(job_id: str):
    """
    Basic UUID format validation to prevent path traversal attacks.
    A valid job_id is a UUID4 string (36 chars, 4 hyphens).
    """
    uuid4_pattern = re.compile(
        r"^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$",
        re.IGNORECASE,
    )
    if not uuid4_pattern.match(job_id):
        raise HTTPException(
            status_code=400,
            detail="Invalid job_id format. Must be a valid UUID4 string.",
        )
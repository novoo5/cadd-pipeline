import asyncio
import shutil
import time
from pathlib import Path
from fastapi import APIRouter
from fastapi.responses import JSONResponse

router = APIRouter()

JOBS_DIR = Path("storage/jobs")

# Track server start time for uptime reporting
_START_TIME = time.time()


@router.get("/health")
async def health_check():
    """
    Primary health check endpoint.
    Hugging Face Spaces pings this to decide if the container is healthy.
    Also used by the frontend to check if the backend is awake before submitting.
    Returns system status including Redis connectivity and disk usage.
    """
    redis_ok   = await _check_redis()
    disk_usage = _get_disk_usage()
    uptime_sec = int(time.time() - _START_TIME)

    status = "healthy" if redis_ok else "degraded"

    return JSONResponse(
        status_code=200 if redis_ok else 503,
        content={
            "status":         status,
            "uptime_seconds": uptime_sec,
            "redis":          "connected" if redis_ok else "unreachable",
            "disk": {
                "used_gb":  disk_usage["used_gb"],
                "free_gb":  disk_usage["free_gb"],
                "total_gb": disk_usage["total_gb"],
            },
            "jobs_on_disk": _count_jobs(),
        },
    )


@router.get("/health/ping")
async def ping():
    """
    Ultra-lightweight keepalive endpoint.
    The frontend can ping this every 4 minutes to prevent HF Space from sleeping.
    (Free HF Spaces sleep after ~15 minutes of inactivity.)
    """
    return {"pong": True, "ts": int(time.time())}


@router.delete("/health/cleanup")
async def cleanup_old_jobs(older_than_hours: int = 24):
    """
    Deletes job directories older than `older_than_hours` hours.
    Run this periodically to free up HF Space storage (5GB free tier limit).
    Can be called from a frontend admin panel or a cron-like ping.
    """
    if not JOBS_DIR.exists():
        return {"deleted": 0, "message": "No jobs directory found."}

    cutoff = time.time() - (older_than_hours * 3600)
    deleted = 0
    errors  = []

    for job_dir in JOBS_DIR.iterdir():
        if not job_dir.is_dir():
            continue
        try:
            mtime = job_dir.stat().st_mtime
            if mtime < cutoff:
                shutil.rmtree(job_dir)
                deleted += 1
        except Exception as e:
            errors.append(str(e))

    return {
        "deleted": deleted,
        "errors":  errors,
        "message": f"Deleted {deleted} job(s) older than {older_than_hours}h.",
    }


@router.get("/health/askcos")
async def check_askcos():
    """
    Checks if the ASKCOS local service is running on port 9100.
    Used by the frontend to warn the user if retrosynthesis is unavailable.
    """
    import httpx
    try:
        async with httpx.AsyncClient(timeout=5.0) as client:
            r = await client.get("http://localhost:9100/api/v2/celery-status/")
        return {
            "askcos": "running" if r.status_code == 200 else "error",
            "http_status": r.status_code,
        }
    except Exception as e:
        return {
            "askcos": "unreachable",
            "error":  str(e),
            "note":   "ASKCOS may still be starting up. Wait 60–90s after container start.",
        }


# ── Internal helpers ─────────────────────────────────────────────────────────

async def _check_redis() -> bool:
    """Pings localhost Redis. Returns True if responsive."""
    try:
        import redis.asyncio as aioredis
        r = aioredis.Redis(host="localhost", port=6379)
        await r.ping()
        await r.aclose()
        return True
    except Exception:
        return False


def _get_disk_usage() -> dict:
    """Returns current disk usage stats for the HF Space storage directory."""
    try:
        usage = shutil.disk_usage("/")
        return {
            "used_gb":  round(usage.used  / 1e9, 2),
            "free_gb":  round(usage.free  / 1e9, 2),
            "total_gb": round(usage.total / 1e9, 2),
        }
    except Exception:
        return {"used_gb": 0, "free_gb": 0, "total_gb": 0}


def _count_jobs() -> int:
    """Counts how many job directories exist on disk."""
    if not JOBS_DIR.exists():
        return 0
    return sum(1 for p in JOBS_DIR.iterdir() if p.is_dir())

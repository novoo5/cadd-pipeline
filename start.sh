#!/bin/bash
set -e

echo "===== Application Startup at $(date -u '+%Y-%m-%d %H:%M:%S') ====="

echo "========================================"
echo "  CADD Pipeline Backend Starting Up"
echo "========================================"

export PYTHONPATH=/workspace

echo "[1/3] Starting Redis server..."
redis-server --daemonize yes
sleep 2
if redis-cli ping > /dev/null 2>&1; then
    echo "      Redis: OK"
else
    echo "      Redis: FAILED TO START"
    exit 1
fi

echo "[2/3] Starting ARQ pipeline worker..."
arq backend.worker.WorkerSettings > /workspace/backend/storage/jobs/arq_worker.log 2>&1 &
echo "      ARQ worker PID: $!"
sleep 2

echo "[3/3] Starting FastAPI on port 7860..."
echo "========================================"
exec uvicorn backend.main:app --host 0.0.0.0 --port 7860

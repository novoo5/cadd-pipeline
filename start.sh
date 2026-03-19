#!/bin/bash
set -e

echo "========================================"
echo "  CADD Pipeline Backend Starting Up"
echo "========================================"

# Let Python know the current directory is the root for imports
export PYTHONPATH=/app:$PYTHONPATH

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
# Run arq using just the module name instead of backend.worker
arq worker.WorkerSettings > storage/jobs/arq_worker.log 2>&1 &
echo "      ARQ worker PID: $!"
sleep 2

echo "[3/3] Starting FastAPI on port 7860..."
echo "========================================"
# Run uvicorn using just main:app instead of backend.main:app
exec uvicorn main:app --host 0.0.0.0 --port 7860

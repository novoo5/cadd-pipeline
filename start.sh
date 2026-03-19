#!/bin/bash
set -e

echo "========================================"
echo "  CADD Pipeline Backend Starting Up"
echo "========================================"

echo "[1/3] Starting Redis server..."
redis-server --daemonize yes --loglevel warning
sleep 1
redis-cli ping > /dev/null 2>&1 && \
    echo "      Redis: OK" || \
    echo "      Redis: FAILED"

echo "[2/3] Starting ARQ pipeline worker..."
python -m arq backend.worker.WorkerSettings &
ARQ_PID=$!
echo "      ARQ worker PID: $ARQ_PID"
sleep 2

echo "[3/3] Starting FastAPI on port 7860..."
echo "========================================"
exec uvicorn backend.main:app \
    --host 0.0.0.0 \
    --port 7860 \
    --workers 1 \
    --log-level info \
    --timeout-keep-alive 75

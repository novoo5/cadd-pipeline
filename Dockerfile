FROM python:3.11-slim

# ── System dependencies ───────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git \
    openbabel libopenbabel-dev \
    libxrender1 libxext6 libx11-6 \
    redis-server \
    build-essential libssl-dev libffi-dev libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# ── GNINA docking binary ──────────────────────────────────────────────────────
RUN wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina \
    -O /usr/local/bin/gnina \
    && chmod +x /usr/local/bin/gnina

# ── PyTorch (CPU) ─────────────────────────────────────────────────────────────
RUN pip install --no-cache-dir torch==2.2.2 \
    --index-url https://download.pytorch.org/whl/cpu

# ── CRITICAL: Pin numpy BEFORE dgl so dgl doesn't pull numpy 2.x ─────────────
RUN pip install --no-cache-dir "numpy==1.26.4"

# ── DGL ───────────────────────────────────────────────────────────────────────
RUN pip install --no-cache-dir dgl==2.1.0 \
    -f https://data.dgl.ai/wheels/repo.html

# ── Copy requirements and install ─────────────────────────────────────────────
# Code lives at /workspace/backend so "from backend.x import y" resolves correctly
WORKDIR /workspace
COPY requirements.txt .
RUN pip install --no-cache-dir --use-deprecated=legacy-resolver -r requirements.txt

# ── AiZynthFinder model download ──────────────────────────────────────────────
RUN python -c "from aizynthfinder.tools.download_public_data import main; main('/root/.aizynthfinder')" || echo "[AiZynth] Model download failed — retrosynthesis will use fallback"

# ── Copy application code into /workspace/backend ────────────────────────────
COPY . /workspace/backend

# ── Storage and startup ───────────────────────────────────────────────────────
RUN mkdir -p /workspace/backend/storage/jobs
RUN chmod +x /workspace/backend/start.sh

ENV PYTHONPATH=/workspace
WORKDIR /workspace

CMD ["/workspace/backend/start.sh"]

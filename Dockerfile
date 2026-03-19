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
RUN pip install --no-cache-dir \
    torch==2.2.2 \
    --index-url https://download.pytorch.org/whl/cpu

# ── CRITICAL: Pin numpy BEFORE dgl so dgl doesn't pull numpy 2.x ─────────────
RUN pip install --no-cache-dir "numpy==1.26.4"

# ── DGL (Graph Neural Networks for binding pre-filter) ────────────────────────
RUN pip install --no-cache-dir \
    dgl==2.1.0 \
    -f https://data.dgl.ai/wheels/repo.html

# ── App dependencies ──────────────────────────────────────────────────────────
WORKDIR /app
COPY requirements.txt .

# We use the legacy resolver to bypass the RDKit/Pandas mutual exclusivity conflict
RUN pip install --no-cache-dir --use-deprecated=legacy-resolver -r requirements.txt

# ── AiZynthFinder model download (optional — has fallback if it fails) ────────
RUN python -c "from aizynthfinder.tools.download_public_data import main; main('/root/.aizynthfinder')" || echo "[AiZynth] Model download failed — retrosynthesis will use fallback"

# ── Copy application code ─────────────────────────────────────────────────────
COPY . .

# ── Storage directory for job files ──────────────────────────────────────────
RUN mkdir -p storage/jobs

# ── Startup script ────────────────────────────────────────────────────────────
RUN chmod +x /app/start.sh
CMD ["/app/start.sh"]

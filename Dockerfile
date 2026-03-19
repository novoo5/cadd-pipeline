FROM python:3.11-slim

LABEL maintainer="cadd-pipeline"
LABEL description="CADD Pipeline: Analogue generation → ADMET → GNINA docking → Retrosynthesis"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# ── System dependencies ───────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    openbabel \
    libopenbabel-dev \
    libxrender1 \
    libxext6 \
    libx11-6 \
    redis-server \
    build-essential \
    libssl-dev \
    libffi-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# ── Install GNINA 1.3 binary ──────────────────────────────────────────────────
RUN wget -q \
    https://github.com/gnina/gnina/releases/download/v1.3/gnina \
    -O /usr/local/bin/gnina \
    && chmod +x /usr/local/bin/gnina

# ── Install PyTorch CPU-only ──────────────────────────────────────────────────
RUN pip install --no-cache-dir \
    torch==2.2.2 \
    --index-url https://download.pytorch.org/whl/cpu

# ── Install DGL CPU-only ──────────────────────────────────────────────────────
RUN pip install --no-cache-dir \
    dgl==2.1.0 \
    -f https://data.dgl.ai/wheels/repo.html

# ── Install Python dependencies ───────────────────────────────────────────────
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ── Download AiZynthFinder pretrained models ──────────────────────────────────
RUN python -c "\
    from aizynthfinder.tools.download_public_data import main; \
    main('/root/.aizynthfinder')" || \
    echo "[AiZynth] Model download failed — retrosynthesis will use fallback"

# ── Copy application code ─────────────────────────────────────────────────────
COPY . .

# ── Create storage directories ────────────────────────────────────────────────
RUN mkdir -p storage/jobs

# ── Copy and set permissions on startup script ───────────────────────────────
RUN chmod +x /app/start.sh

# ── HF Spaces requires port 7860 ─────────────────────────────────────────────
EXPOSE 7860

# ── Health check ─────────────────────────────────────────────────────────────
HEALTHCHECK --interval=30s --timeout=10s --start-period=90s --retries=3 \
    CMD curl -f http://localhost:7860/api/v1/health/ping || exit 1

CMD ["/app/start.sh"]

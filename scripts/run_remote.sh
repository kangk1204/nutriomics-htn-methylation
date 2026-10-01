#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export JOBLIB_TEMP_FOLDER="$PWD/cache/joblib"
mkdir -p "$JOBLIB_TEMP_FOLDER"
METHYLATION_PYTHON="${METHYLATION_PYTHON:-/home/keunsoo/projects/nutriomics-2026-final/.venv/bin/python}"
"$METHYLATION_PYTHON" -m pytest -q
"$METHYLATION_PYTHON" -m nutriomics_methylation.cli train --cores 12 --spaces common full --methods nested5 loocv
"$METHYLATION_PYTHON" -m nutriomics_methylation.cli evaluate --bootstrap-replicates 1000
"$METHYLATION_PYTHON" -m nutriomics_methylation.cli export

#!/bin/bash
#SBATCH --job-name=eps_shock
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=eps_shock_%A_%a.out
#SBATCH --error=eps_shock_%A_%a.err

## Submit from the directory containing:
##   shock_on_eps.R
##   Data_Setting.Rdata
##
## Command:
##   sbatch shock_on_eps.sbatch

set -euo pipefail

mkdir -p results_eps_shock

## Clusters can differ in module names. If your cluster uses a versioned
## R module, replace "module load R" with that module name.
if command -v module >/dev/null 2>&1; then
  module purge || true
  module load R || true
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript was not found. Please load the correct R module in this sbatch file."
  exit 1
fi

export INPUT_RDATA="Data_Setting.Rdata"
export OUTPUT_DIR="results_eps_shock"
export RUN="10000"
export LAG="20"
export BASE_SEED="1130"

export USE_MARGINAL_SD="false"

Rscript shock_on_eps.R "${SLURM_ARRAY_TASK_ID}" "${RUN}" "${LAG}"

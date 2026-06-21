#!/bin/bash
#SBATCH --job-name=0621_HAC
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=0621_HAC_%A_%a.out
#SBATCH --error=0621_HAC_%A_%a.err

## Submit from the directory containing:
##   0621_HAC_satterthwaite_prewhite_cluster.R
##   1108_PermSimu_Setting_c.Rdata
##
## Command:
##   sbatch 0621_HAC_array.sbatch

set -euo pipefail

mkdir -p results_0621

## Yale clusters can differ in module names. If your cluster uses a versioned
## R module, replace "module load R" with that module name.
if command -v module >/dev/null 2>&1; then
  module purge || true
  module load R || true
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript was not found. Please load the correct R module in this sbatch file."
  exit 1
fi

export INPUT_RDATA="1108_PermSimu_Setting_c.Rdata"
export OUTPUT_DIR="results_0621"
export RUN="10000"
export LAG="20"
export BASE_SEED="1130"

## Keep FALSE to reproduce the convention in your 0621.R.
## Set TRUE only if you want sd(eps_perp_ori_wide) and sd(Xperp_ori_wide)
## to be treated as marginal AR(1) standard deviations rather than innovation SDs.
export USE_MARGINAL_SD="false"

Rscript 0621_HAC_satterthwaite_prewhite_cluster.R "${SLURM_ARRAY_TASK_ID}" "${RUN}" "${LAG}"

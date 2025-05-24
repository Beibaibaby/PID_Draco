#!/bin/bash
#SBATCH --job-name=FIT6
#SBATCH --account=pi-bdoiron
#SBATCH --partition=caslake
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16           # each task gets 16 cores
#SBATCH --mem=64G
#SBATCH --array=0-5                  # 6 combos (see list below)
#SBATCH --output=FIT_%A_%a.txt    # <-- stdout
#SBATCH --error=FIT_%A_%a.err    # <-- stderr

module load matlab

# -------- pair/label lookup table ----------
declare -a COMBOS=(
    "LIP FEF direction"
    "LIP FEF category"
    "FEF SC  direction"
    "FEF SC  category"
    "LIP SC  direction"
    "LIP SC  category"
)
set -- ${COMBOS[$SLURM_ARRAY_TASK_ID]}
SRC=$1 ; DST=$2 ; LABEL=$3
# -------------------------------------------

echo "[$(date)]  Task $SLURM_ARRAY_TASK_ID ➜ $SRC->$DST ($LABEL)"

cd /project/bdoiron/dracoxu/PID_Draco/fit_3pair_pipeline
matlab -nodisplay -r "\
  addpath(genpath('/project/bdoiron/draco/MINT')); \
  addpath(genpath(pwd)); \
  run_fit_combo('$SRC','$DST','$LABEL'); \
  matlab -nodisplay -r "run_postPlots";\
  quit"

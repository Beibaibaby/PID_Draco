#!/bin/bash
#SBATCH --job-name=FIT-tau
#SBATCH --account=pi-bdoiron
#SBATCH --partition=caslake
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-29                          # 6 pairs × 5 τ = 30
#SBATCH --output=logs/FIT_%A_%a.txt    
#SBATCH --error=logs/FIT_%A_%a.err

module load matlab

# ------------------------------------------------------------------
# create the sub-directory for this array task (%a) if it doesn’t exist
mkdir -p "logs/tau${SLURM_ARRAY_TASK_ID}"

# your task-lookup table etc. unchanged …
declare -a COMBOS=(
    "LIP FEF direction"
    "LIP FEF category"
    "FEF SC  direction"
    "FEF SC  category"
    "LIP SC  direction"
    "LIP SC  category"
)
set -- ${COMBOS[$SLURM_ARRAY_TASK_ID]}
SRC=$1 ; DST=$2 ; LABEL=$3 ; TAU=$(( ($SLURM_ARRAY_TASK_ID / 6) + 1 ))  # 1..5

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID → $SRC->$DST, $LABEL, τ=$TAU"

cd /project/bdoiron/dracoxu/PID_Draco/fit_3pair_pipeline

matlab -nodisplay -r "\
 addpath(genpath('/project/bdoiron/draco/MINT')); \
 addpath(genpath(pwd)); \
 run_fit_combo('$SRC','$DST','$LABEL',$TAU); \
 quit"

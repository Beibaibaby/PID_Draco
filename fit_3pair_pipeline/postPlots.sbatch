#!/bin/bash
#SBATCH --job-name=FIT-plot
#SBATCH --account=pi-bdoiron
#SBATCH --partition=caslake
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/postPlots_%j.txt
#SBATCH --error =logs/postPlots_%j.err

module load matlab
cd /project/bdoiron/dracoxu/PID_Draco/fit_3pair_pipeline
matlab -nodisplay -r "run_postPlots"

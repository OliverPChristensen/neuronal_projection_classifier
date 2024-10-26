#!/bin/bash

#SBATCH --job-name=spatial_custering_analysis
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/running_sbatch_logs/%j.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/running_sbatch_logs/%j_error.log

#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

Rscript scripts/spatial_custering_analysis.R -o $1 -p $2 -sop $3 -rps $4 -rdp $5

pid=$!
wait $pid
echo Moving log file and error log file
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/running_sbatch_logs/"$SLURM_JOB_ID".log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/spatial_custering_analysis.log
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/running_sbatch_logs/"$SLURM_JOB_ID"_error.log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/spatial_custering_analysis_error.log
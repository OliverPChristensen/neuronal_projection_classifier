#!/bin/bash

#SBATCH --job-name=generate_spatial_object
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_error.log

#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

Rscript scripts/spatial_projection_call.R -o $1 -sop $2 -tp $3 -cp $4 -ip $5 -sp $6 -ri $7 -si $8

pid=$!
wait $pid
echo Moving log file and error file
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID".log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/"$7"/"$8"_spatial_projection_call.log
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID"_error.log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/"$7"/"$8"_spatial_projection_call_error.log
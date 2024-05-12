#!/bin/bash

#SBATCH --job-name=generate_spatial_object
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_error.log

#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

Rscript scripts/generate_spatial_object.R -o $1 -tp $2 -cp $3 -ip $4 -rp $5 -ri $6 -si $7

pid=$!
wait $pid
echo Moving logfile
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID".log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/"$6"/"$7"_generate_spatial_object.log
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID"_error.log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/"$6"/"$7"_generate_spatial_object_error.log
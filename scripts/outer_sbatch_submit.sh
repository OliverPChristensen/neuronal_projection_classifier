#!/bin/bash

#SBATCH --job-name=segmentation
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/segmentation_outer%a.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/segmentation_outer_error%a.log
#SBATCH --array=1-3

##SBATCH --cpus-per-task=1
##SBATCH --mem=16G
#SBATCH --time=20-00:00:00

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier

run_indices=("spatial1" "spatial2" "spatial3")
run_index=${run_indices[$SLURM_ARRAY_TASK_ID-1]}

sbatch scripts/inner_sbatch_submit.sh $SLURM_ARRAY_TASK_ID $run_index

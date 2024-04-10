#!/bin/bash

#SBATCH --job-name=segmentation
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner8.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner_error8.log

#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=20-00:00:00

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

python scripts/prepare_data_for_seurat.py -ri spatial1 -si D2
#!/bin/bash

#SBATCH --job-name=segmentation
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner%a.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner_error%a.log
#SBATCH --array=1-8

#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=20-00:00:00

section_indices=("A1" "B1" "C1" "D1" "A2" "B2" "C2" "D2")
polyT_section_indices=("W0" "W1" "W2" "W3" "W4" "W5" "W6" "W7")
section_index=${section_indices[$SLURM_ARRAY_TASK_ID-1]}
polyT_section_index=${polyT_section_indices[$SLURM_ARRAY_TASK_ID-1]}

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

python scripts/prepare_data_for_seurat.py -ri $2 -si $section_index -psi $polyT_section_index -dfm DAPI -pfm Cy5
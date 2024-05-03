#!/bin/bash

#SBATCH --job-name=segmentation
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner%a.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner_error%a.log
#SBATCH --array=1-8

#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=20-00:00:00

section_indices=("A1" "B1" "C1" "D1" "A2" "B2" "C2" "D2")
indices_converter=("W0" "W1" "W2" "W3" "W4" "W5" "W6" "W7")
section_index=${section_indices[$SLURM_ARRAY_TASK_ID-1]}
index_converter=${indices_converter[$SLURM_ARRAY_TASK_ID-1]}

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

python scripts/cell_segmentation_pipeline.py -t combined -ri $2 -si $section_index -ic2 $index_converter -fm DAPI -fm2 Cy5 -nt cell_seg

pid=$!
wait $pid
echo Moving logfile
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner%a.log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/"$2"/"$section_index"_cell_seg.log
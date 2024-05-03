#!/bin/bash

#SBATCH --job-name=segmentation
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner8.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_segmentation_inner_error8.log

#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=20-00:00:00

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier
echo Xenium_test_resolve
python scripts/cell_segmentation_pipeline.py -t single -ri spatial3 -si B1 -fm DAPI -nt xenium_test_resolve

pid=$!
wait $pid
echo Moving logfile
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID"_segmentation_inner8.log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/"$2"/"$section_index"_cell_seg_xenium_test_resolve.log


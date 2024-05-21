#!/bin/bash

#SBATCH --job-name=presentation
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_error.log

#SBATCH --cpus-per-task=1
#SBATCH --mem=75G
#SBATCH --time=20-00:00:00

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

python scripts/cell_segmentation_pipeline.py -o presentations/shared_lab_meeting_24052024/segmentations -sp raw_data/spatial1/Spatial1_serie2_serie2_A1-1_DAPI.tiff -sp2 raw_data/spatial1/Panorama_Spatial1_serie2_W0A1_Cy5-class_R11_.tiff -ri spatial1 -si A1 -nt combined -gf False

pid=$!
wait $pid
echo Moving logfile
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID".log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/presentations/shared_lab_meeting_24052024/segmentations/segmentation_combined_no_gridfill.log
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID"_error.log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/presentations/shared_lab_meeting_24052024/segmentations/segmentation_combined_no_gridfill_error.log

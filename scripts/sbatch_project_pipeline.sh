#!/bin/bash

#SBATCH --job-name=project_pipeline
#SBATCH --output=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j.log
#SBATCH --error=/projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/%j_error.log

#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=20-00:00:00

cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

snakemake -s scripts/project_pipeline.snakefile -j 30 --latency-wait 60 -R merge_qc_norm_integrate_spatial_object

pid=$!
wait $pid
echo Moving logfile
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID".log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/project_pipeline.log
mv /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/processed_data/segmentation_logs/"$SLURM_JOB_ID"_error.log /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier/project_pipeline_error.log
cd /projects/cbmr_shared/people/wqc597/neuronal_projection_classifier
source activate_project_neuronal_projection_classifier

sbatch scripts/outer_sbatch_submit.sh

Rscript scripts/generate_spatial_objects.R -ri spatial1 -fmt results -fmr rois
Rscript scripts/generate_spatial_objects.R -ri spatial2 -fmt results -fmr rois
Rscript scripts/generate_spatial_objects.R -ri spatial3 -fmt results -fmr rois
# Snakefile
import os
run_indices=("spatial1", "spatial2", "spatial3")
section_indices=("A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2")
indices_converter=("W0", "W1", "W2", "W3", "W4", "W5", "W6", "W7")

# List of input files with full paths
def fetch_file_path(path_to_staining_folder, file_match1,file_match2):
    files = os.listdir(path_to_staining_folder)
    matching_files = [path_to_staining_folder + "/" + file for file in files if file_match1 in file and file_match2 in file]
    return matching_files

ic = dict(zip(section_indices, indices_converter))

rule all:
    input:
        expand("processed_data/{run_index}/spatial_object_merged.qs", run_index=run_indices),
        #expand("processed_data/{run_index}/{section_index}_spatial_object_projection_call.qs", run_index = run_indices, section_index = section_indices),
        expand("processed_data/{run_index}/{section_index}_cell_seg_imagej_rois.zip", run_index = run_indices, section_index = section_indices),
        expand("processed_data/{run_index}/{section_index}_sun1_imagej_rois.zip", run_index = run_indices, section_index = section_indices),
        expand("processed_data/{run_index}/{section_index}_sun1_area_correction_rois.npz", run_index = run_indices, section_index = section_indices)

rule spatial_merge:
    input:
        spatial_objects = expand("processed_data/{run_index}/{section_index}_spatial_object_projection_call.qs",run_index = run_indices, section_index=section_indices)
        #spatial_objects = processed_data/{run_index}/{section_index}_spatial_object_projection_call.qs"
    output:
        "processed_data/{run_index}/spatial_object_merged.qs"
    shell:
        "Rscript scripts/spatial_merge.R -o processed_data/{run_index} -f {input.spatial_objects}"

rule spatial_projection_call:
    input: 
        spatial_object = "processed_data/{run_index}/{section_index}_spatial_object_raw.qs",
        transcripts = lambda wildcards: fetch_file_path("raw_data/" + wildcards.run_index, "results",wildcards.section_index),
        cells = "processed_data/{run_index}/{section_index}_cell_seg_cells_rois.npz",
        image = "processed_data/{run_index}/{section_index}_cell_seg_area_correction_rois.npz",
        sun1 = "processed_data/{run_index}/{section_index}_sun1_cells_rois.npz",
    output: "processed_data/{run_index}/{section_index}_spatial_object_projection_call.qs"
    shell:
        "sbatch --wait scripts/sbatch_spatial_projection_call.sh processed_data/{wildcards.run_index} {input.spatial_object} {input.transcripts} {input.cells} {input.image} {input.sun1} {wildcards.run_index} {wildcards.section_index}"

rule generate_spatial_object:
    input: 
        transcripts = lambda wildcards: fetch_file_path("raw_data/" + wildcards.run_index, "results",wildcards.section_index),
        cells = "processed_data/{run_index}/{section_index}_cell_seg_cells_rois.npz",
        image = "processed_data/{run_index}/{section_index}_cell_seg_area_correction_rois.npz",
        regions = "raw_data/{run_index}/{run_index}_{section_index}_regions"
    output: "processed_data/{run_index}/{section_index}_spatial_object_raw.qs"
    shell:
        "sbatch --wait scripts/sbatch_generate_spatial_object.sh processed_data/{wildcards.run_index} {input.transcripts} {input.cells} {input.image} {input.regions} {wildcards.run_index} {wildcards.section_index}"

rule cell_segmentation_pipeline:
    input:
        dapi = lambda wildcards: fetch_file_path("raw_data/" + wildcards.run_index,"DAPI",wildcards.section_index),
        polyt = lambda wildcards: fetch_file_path("raw_data/" + wildcards.run_index,"Cy5",ic[wildcards.section_index])
    output:
        "processed_data/{run_index}/{section_index}_cell_seg_cells_rois.npz",
        "processed_data/{run_index}/{section_index}_cell_seg_imagej_rois.zip",
        "processed_data/{run_index}/{section_index}_cell_seg_area_correction_rois.npz"
    shell:
        """
        sbatch --wait scripts/sbatch_cell_segmentation_pipeline.sh processed_data/{wildcards.run_index} {input.dapi} {input.polyt} {wildcards.run_index} {wildcards.section_index} cell_seg
        """

rule sun1_segmentation_pipeline:
    input:
        sun1 = lambda wildcards: fetch_file_path("raw_data/" + wildcards.run_index,"GFP",ic[wildcards.section_index])
    output:
        "processed_data/{run_index}/{section_index}_sun1_cells_rois.npz",
        "processed_data/{run_index}/{section_index}_sun1_imagej_rois.zip",
        "processed_data/{run_index}/{section_index}_sun1_area_correction_rois.npz"
    shell:
        "sbatch --wait scripts/sbatch_sun1_segmentation_pipeline.sh processed_data/{wildcards.run_index} {input.sun1} {wildcards.run_index} {wildcards.section_index} sun1"
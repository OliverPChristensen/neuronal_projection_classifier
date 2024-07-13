# Snakefile
import os
import itertools

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
        "spatial_hypothalamus_app/data/shiny_data.qs",
        "plots/spatial_cluster_analysis/chi2_p_values_plot.png",
        "plots/region_projection_analysis/region_projection_plot_bayesian.png",
        "plots/validate_projection_call/check_background_intensity_plot.png"

rule create_shiny_data:
    input:
        spatial_object = "processed_data/spatial_object_preprocessed.qs",
        transcripts = [file for run_index, section_index in itertools.product(run_indices, section_indices) for file in fetch_file_path(f"raw_data/{run_index}", "results", section_index)],
        regions = expand("raw_data/{run_index}/{run_index}_{section_index}_regions", run_index = run_indices, section_index = section_indices),
        bregma_tsv = "processed_data/ROI_Bregma_spatial.tsv"
    output:
        "spatial_hypothalamus_app/data/shiny_data.qs"
    shell:
        """
        sbatch --wait scripts/sbatch_create_shiny_data.sh spatial_hypothalamus_app/data {input.spatial_object} "{input.transcripts}" "{input.regions}" {input.bregma_tsv}
        """

rule spatial_cluster_analysis:
    input:
        spatial_object = "processed_data/spatial_object_preprocessed.qs",
        regions = expand("raw_data/{run_index}/{run_index}_{section_index}_regions", run_index = run_indices, section_index = section_indices)
    output:
        "plots/spatial_cluster_analysis/chi2_p_values_plot.png"
    shell:
        """
        sbatch --wait scripts/sbatch_spatial_clustering_analysis.sh processed_data plots/spatial_cluster_analysis {input.spatial_object} "{input.regions}"
        """

rule region_projection_analysis:
    input:
        spatial_object = "processed_data/spatial_object_preprocessed.qs"
    output:
        "plots/region_projection_analysis/region_projection_plot_bayesian.png"
    shell:
        "sbatch --wait scripts/sbatch_region_projection_analysis.sh processed_data plots/region_projection_analysis {input.spatial_object}"

rule validate_projection_call:
    input:
        spatial_object = "processed_data/spatial_object_preprocessed.qs"
    output:
        "plots/validate_projection_call/check_background_intensity_plot.png"
    shell:
        "sbatch --wait scripts/sbatch_validate_projection_call.sh processed_data plots/validate_projection_call {input.spatial_object}"

rule final_preprocessing:
    input:
        spatial_objects = expand("processed_data/{run_index}/{section_index}_spatial_object_projection_call.qs",run_index = run_indices, section_index=section_indices)
    output:
        "processed_data/spatial_object_preprocessed.qs"
    shell:
        """
        sbatch --wait scripts/sbatch_final_preprocessing.sh processed_data plots/final_preprocessing "{input.spatial_objects}"
        """

rule spatial_projection_call:
    input: 
        spatial_object = "processed_data/{run_index}/{section_index}_spatial_object_raw.qs",
        transcripts = lambda wildcards: fetch_file_path("raw_data/" + wildcards.run_index, "results",wildcards.section_index),
        cells = "processed_data/{run_index}/{section_index}_cell_seg_cells_rois.npz",
        sun1 = "processed_data/{run_index}/{section_index}_sun1_cells_rois.npz",
    output: "processed_data/{run_index}/{section_index}_spatial_object_projection_call.qs"
    shell:
        "sbatch --wait scripts/sbatch_spatial_projection_call.sh processed_data/{wildcards.run_index} {input.spatial_object} {input.transcripts} {input.cells} {input.sun1} {wildcards.run_index} {wildcards.section_index}"

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
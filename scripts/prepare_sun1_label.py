import argparse
import sys
import os
import numpy as np

from skimage import io

sys.path.append(os.getcwd())

print(f"Current working directory: {os.getcwd()}")

from scripts.cell_segmentation_pipeline import current_time, fill_grids, staining_to_rois, save_rois_as_imagej



def load_sun1_staining(path_to_staining_folder,sun1_file_match,section_index,section_index_pairing):
    #Load files in folder with stanings
    files = os.listdir(path_to_staining_folder)
 
    #Fetch all dapi and polyT file names
    sun1_list = np.array([element for element in files if sun1_file_match in element])

    #Find the dapi file name that matches the section index
    sun1_file = sun1_list[np.core.defchararray.find(sun1_list,section_index_pairing[section_index]) > 0].item()

    #Define path to DAPI file and polyT file
    path_to_sun1_file = f"{path_to_staining_folder}/{sun1_file}"
    
    #Load paired dapi and polyT stainings
    sun1_image = io.imread(path_to_sun1_file)

    return sun1_image

def create_parser():
    '''
    Input: NA

    Creates and outputs parser with arguments that allows interaction via CLI

    Output: Parser
    '''

    #Define parser
    parser = argparse.ArgumentParser(description="Segmentation Pipeline")

    #Add argumnent to parser
    parser.add_argument("-fp","--folder-path", type=str, default = 'default', help="Relative path to folder with stainings")
    parser.add_argument("-sfm","--sun1-file-match", type=str, default = '', help="String to be used to specify Sun1 files in stainings folder")
    parser.add_argument("-o","--output-folder-path", type=str, default = './processed_data', help="Relative path to folder to save output")
    parser.add_argument("-si","--section-indices", nargs="+", type=str, default=["A1","B1","C1","D1","A2","B2","C2","D2"], help="List of section indices to segment")
    parser.add_argument("-ic","--index_converter", nargs="+", type=str, default='default', help="List of section indices to use for naming output")
    parser.add_argument("-ri","--run-index", type=str, help="List of section indices to segment")

    return parser

def main():
    '''
    Input: NA

    Main function that defines parser arguments, creates section index pairing dictionary, loades DAPI and polyT image, 
    combines them, runs cellpose on the combined image and creates ROIs from the masks outputted, save the ROIs and generate a quality report of the ROIs

    Output: NA
    '''

    #Create parser and parse arguments
    parser = create_parser()
    args = parser.parse_args()
    
    #Define variables from parser
    sun1_file_match = args.sun1_file_match
    path_to_output_folder = args.output_folder_path
    section_indices = args.section_indices
    index_converter = args.index_converter
    run_index = args.run_index

    if args.folder_path == 'default':
        path_to_staining_folder = f"./raw_data/{run_index}"
    else:
        path_to_staining_folder = args.folder_path

    #Dictionary to convert between naming schemes for different section index depending on polyT spatial indices
    if index_converter == 'default':
        section_index_pairing = dict(zip(section_indices, section_indices))
    else:
        section_index_pairing = dict(zip(section_indices, index_converter))

    #Create folder to output ROIs
    os.makedirs(f'{path_to_output_folder}/{run_index}',exist_ok=True)

    #Loop over spatial indices 
    for i, section_index in enumerate(section_indices):

        #Load sun1 image
        print(f"[{current_time()}] Loading section {section_index} from {path_to_staining_folder} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        image = load_sun1_staining(path_to_staining_folder,sun1_file_match,section_index,section_index_pairing)

        #Fill grids in stainings
        print(f"[{current_time()}] Filling grids on section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        image = fill_grids(image)

        #Input the combined image to staining_to_rois, which returns cell ROIs for the image
        print(f"[{current_time()}] Running cellpose on section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        rois = staining_to_rois(image)

        #Save the ROIs
        print(f"[{current_time()}] Saving Sun1 rois for section {section_index} to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        np.savez(f'{path_to_output_folder}/{run_index}/{section_index}_sun1_rois.npz', **rois)

        #Save the ROIs in ImageJ format
        print(f"[{current_time()}] Saving Sun1 rois for section {section_index} in ImageJ format to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        save_rois_as_imagej(rois, path_to_output_folder, run_index, section_index, "sun1")

if __name__ == "__main__":
    main()



import os
import sys
from datetime import datetime
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from skimage import io, measure

from cellpose import io
from cellpose import models

from roifile import ImagejRoi, roiwrite

def current_time():
    '''
    Allows the current time to be added to all prints to the log, when script is run through sbatch
    Output: Current time
    '''
    cur_time = datetime.now()
    return cur_time.strftime("%Y-%m-%d %H:%M:%S")

def staining_to_rois(image):
    '''
    Input: Numpy array of staining image

    Runs cellpose on the combined staining image
    Converts masks to dictionary of ROIs

    Output: Dictionary of ROIs for the segemnetation of the stainings
    '''

    #Define cellpose model
    model = models.Cellpose(gpu=True, model_type = 'cyto')

    #Print device
    print(f"[{current_time()}]: Device is {model.device}")
    sys.stdout.flush()

    #Generate masks for image
    masks, _, _, diams = model.eval(image, diameter=None, channels=[0,0],
                                         flow_threshold=0.4, do_3D=False)
    
    #Convert masks from image into ROIs and put them into a dictionary
    rois = {}
    for label in range(1, np.max(masks) + 1):
        mask = (masks == label).astype(np.uint8)
        contour = measure.find_contours(mask, 0.5)
        rois[f"Cell{label}"] = contour[0][:,[1,0]]

    return rois

def quality_report(image,rois,run,section):
    '''
    Input: Numpy array of staining image, dictionary of ROIs, spatial run index, section ID

    Creates and saves multiple plots of the staining image with ROIs overload. 
    The image is split into smaller images to be able to open them with sufficient resolution in VS Code

    Output: NA
    '''

    #Define plot and add image to plot
    fig, ax = plt.subplots()
    _ = ax.imshow(image)
    
    #Add ROIs to plot
    for roi in rois:
        patch = patches.Polygon(rois[roi], closed=True, edgecolor='r', facecolor='none')
        _ = ax.add_patch(patch)

    #Define split on each dimension of the image
    image_split_num = 50

    #Calculate width of each split
    xdim = int(image.shape[1]/image_split_num)
    ydim = int(image.shape[0]/image_split_num)
    
    #Create folder to save plots
    os.makedirs(f"./plots/segmentation_report/{run}/{section}",exist_ok=True)

    #Generate all plots as specified by the splits. Empty images are ignored
    for i in range(image_split_num):
        for j in range(image_split_num):
            
            #Checks if cropped image is empty. If not, then save plot
            if np.max(image[ydim*j:ydim*(j + 1),xdim*i:xdim*(i + 1)]) > 0:
                _ = ax.set_xlim(xdim*i, xdim*(i + 1))
                _ = ax.set_ylim(ydim*j, ydim*(j + 1))

                plt.savefig(f"./plots/segmentation_report/{run}/{section}/{section}_quality_report_{i}-{j}.png")

def save_rois_as_imagej(rois, path_to_output_folder, run_index, section_index):
    '''
    Input: ROIs in a dictionary, path to folder input, run index, section index

    Converts the ROIs into imageJ format and saves them in a zip file

    Output: NA
    '''

    #Convert ROIs into imageJ format
    rois_imJ = [ImagejRoi.frompoints(outline) for _, outline in rois.items()]

    #Create filename for the ROIs
    file_name = f'{path_to_output_folder}/{run_index}/{section_index}_imagej_rois.zip'

    #Check if the file already exists. This is necessary because roiwrite will append rois to an already existing file instead of overwriting it
    if os.path.exists(file_name):
            os.remove(file_name)

    #Save the ROIs as a zip file
    roiwrite(file_name, rois_imJ)

def load_dapi_polyT_pair(path_to_staining_folder,dapi_file_match,polyT_file_match,section_index,section_index_pairing):
    '''
    Input: Path to staining folder, DAPI file match string, polyT file match string, section index, dictionary with section index pairing

    Loads DAPI path according to section index, finds corresponding polyT path using section index pairing, loads the DAPI and polyT image

    Output: DAPI image as numpy array, polyT image as numpy array
    '''
    #Load files in folder with stanings
    files = os.listdir(path_to_staining_folder)

    #Fetch all dapi and polyT file names
    dapi_list = np.array([element for element in files if dapi_file_match in element])
    polyT_list = np.array([element for element in files if polyT_file_match in element])

    #Find the dapi file name that matches the section index
    dapi_file = dapi_list[np.core.defchararray.find(dapi_list,section_index) > 0].item()

    #Fetch the alternative section index paired to the section index
    alt_section_index = section_index_pairing[section_index]

    #Find the polyT file name that matches the alternative section index
    polyT_file = polyT_list[np.core.defchararray.find(polyT_list, alt_section_index) > 0].item()    

    #Define path to DAPI file and polyT file
    path_to_dapi_file = f"{path_to_staining_folder}/{dapi_file}"
    path_to_polyT_file = f"{path_to_staining_folder}/{polyT_file}"
    
    #Load paired dapi and polyT stainings
    dapi_image = io.imread(path_to_dapi_file)
    polyT_image = io.imread(path_to_polyT_file)

    return dapi_image, polyT_image

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
    parser.add_argument("-dfm","--dapi-file-match", type=str, default = '', help="String to be used to specify DAPI files in stainings folder")
    parser.add_argument("-pfm","--polyT-file-match", type=str, default = '', help="String to be used to specify polyT files in stainins folder")
    parser.add_argument("-psi","--polyT-spatial-indices", nargs="+", type=str, default = 'default', help="String to be used to specify polyT files in stainins folder")
    parser.add_argument("-o","--output-folder-path", type=str, default = './processed_data', help="Relative path to folder to save output")
    parser.add_argument("-si","--section-indices", nargs="+", type=str, default=["A1","B1","C1","D1","A2","B2","C2","D2"], help="List of section indices to segment")
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
    dapi_file_match = args.dapi_file_match
    polyT_file_match = args.polyT_file_match
    path_to_output_folder = args.output_folder_path
    section_indices = args.section_indices
    polyT_spatial_indices = args.polyT_spatial_indices
    run_index = args.run_index

    if args.folder_path == 'default':
        path_to_staining_folder = f"./raw_data/{run_index}"
    else:
        path_to_staining_folder = args.folder_path

    #Dictionary to convert between naming schemes for different section index depending on polyT spatial indices
    if polyT_spatial_indices == 'default':
        section_index_pairing = dict(zip(section_indices, section_indices))
    else:
        section_index_pairing = dict(zip(section_indices, polyT_spatial_indices))

    #Loop over spatial indices 
    for i, section_index in enumerate(section_indices):

        #Load DAPI image and polyT image according to section index
        print(f"[{current_time()}]: Loading section {section_index} from {path_to_staining_folder} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        dapi_image, polyT_image = load_dapi_polyT_pair(path_to_staining_folder,dapi_file_match,polyT_file_match,section_index,section_index_pairing)

        #Combine dapi and polyT stainings by using average pixel value
        image = (dapi_image + polyT_image)/2

        #Input the combined image to staining_to_rois, which returns cell ROIs for the image
        print(f"[{current_time()}]: Running cellpose on section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        rois = staining_to_rois(image)

        #Save the ROIs
        print(f"[{current_time()}]: Saving rois for section {section_index} to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        np.savez(f'{path_to_output_folder}/{run_index}/{section_index}_rois.npz', **rois)

        #Save the ROIs in ImageJ format
        print(f"[{current_time()}]: Saving rois for section {section_index} in ImageJ format to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        save_rois_as_imagej(rois, path_to_output_folder, run_index, section_index)
        
        """
        #Generate quility report of the segmentation using quality_report
        print(f"[{current_time()}]: Generating quality report for section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        quality_report(image, rois, run_index, section_index)
        """

if __name__ == "__main__":
    main()
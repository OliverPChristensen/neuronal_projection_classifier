import os
import sys
from datetime import datetime
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from skimage import io, measure

#from cellpose import io
from cellpose import models

from roifile import ImagejRoi, roiwrite

import cv2

def current_time():
    '''
    Allows the current time to be added to all prints to the log
    Output: Current time
    '''
    cur_time = datetime.now()
    return cur_time.strftime("%Y-%m-%d %H:%M:%S")


def load_staining(path_to_staining_folder,file_match,section_index,section_index_pairing):
    '''
    Input: Path to staining folder, File match string, Section index, dictionary with section index pairing

    Loads image according to file match and section index

    Output: Image as numpy array
    '''
    #Load files in folder with stanings
    files = os.listdir(path_to_staining_folder)
    
    #Fetch all file names that match the file match
    matching_files = np.array([file for file in files if file_match in file])

    #Find the file name that matches the section index
    file_name = matching_files[np.core.defchararray.find(matching_files,section_index_pairing[section_index]) > -1].item()

    #Define path to staining file
    path_to_file = f"{path_to_staining_folder}/{file_name}"
    
    #Load staining file
    image = io.imread(path_to_file)

    return image


def generate_rois_for_area_correction(image):
    '''
    Input: Image before grid fill

    Generates ROIs for all non zero areas in the image before grid fill. 
    This can be used to correct for overestimation of areas as a result of grid fill

    Output: Dictionary with ROIs of non zero pixel value areas in the image
    '''

    #Generate mask(s) with all non-zero pixel values of the image
    mask = (image != 0).astype(np.uint8)

    #Find ROIs for the mask(s)
    contour = measure.find_contours(mask, 0.5)

    #Convert the ROIs to a dictionary format and invert axis of coordinates
    rois = {}
    for i, cont in enumerate(contour):
        rois[str(i)] = cont[:,[1,0]]

    return rois

def fill_grids(img_array, box_size = 5, nloops = 100):
    '''
    Input: 1 channel image as numpy array, size of kernel to use, number of rounds to apply gaussian blur

    MODIFIED VERSION OF MINDAGAP: https://github.com/ViriatoII/MindaGap/tree/main 
    Fills grid locations (pixel val 0) with values from neighbourhood through a gaussian kernel 

    Output: Grid filled image as numpy array
    '''

    # Grid coordinates and a copy of image
    grid_coords = img_array == 0
    im_copy = img_array.copy()   
   
    # Make grid pixels have the minimum value of image (excluding 0 of grid)
    im_copy[grid_coords] = min(img_array[img_array > 0].flatten()) 
   
    # Create a blurred image and replace original grid positions by new blur value. Iterate
    for i in range(nloops):

        blur_img = cv2.GaussianBlur(im_copy,(box_size,box_size), 0)   # Gaussian kernel
        im_copy[grid_coords] = blur_img[grid_coords]

    return(im_copy)


def staining_to_rois(image,diameter):
    '''
    Input: Numpy array of staining image, diameter to use for segmentation

    Runs cellpose on the combined staining image
    Converts masks to dictionary of ROIs

    Output: Dictionary of ROIs for the segemnetation of the stainings
    '''

    #Define cellpose model
    model = models.Cellpose(gpu=True, model_type = 'cyto')

    #Print device
    print(f"[{current_time()}] Device is {model.device}")
    sys.stdout.flush()

    #Generate masks for image
    masks, _, _, diams = model.eval(image, diameter=diameter, channels=[0,0],
                                         flow_threshold=0.4, do_3D=False)
    
    print(f"[{current_time()}] Diameter used for the segmentation: {diams}")
    sys.stdout.flush()
    
    #Convert masks from image into ROIs and put them into a dictionary
    rois = {}
    for label in range(1, np.max(masks) + 1):
        mask = (masks == label).astype(np.uint8)
        contour = measure.find_contours(mask, 0.5)
        rois[f"Cell{label}"] = contour[0][:,[1,0]]

    return rois


def save_rois_as_imagej(rois, path_to_output_folder, run_index, section_index, name_tag):
    '''
    Input: ROIs in a dictionary, path to folder input, run index, section index

    Converts the ROIs into imageJ format and saves them in a zip file

    Output: NA
    '''

    #Convert ROIs into imageJ format
    rois_imJ = [ImagejRoi.frompoints(outline) for _, outline in rois.items()]

    #Create filename for the ROIs
    file_name = f'{path_to_output_folder}/{run_index}/{section_index}_{name_tag}_imagej_rois.zip'

    #Check if the file already exists. This is necessary because roiwrite will append rois to an already existing file instead of overwriting it
    if os.path.exists(file_name):
            os.remove(file_name)

    #Save the ROIs as a zip file
    roiwrite(file_name, rois_imJ)


def create_parser():
    '''
    Input: NA

    Creates and outputs parser with arguments that allows interaction via CLI

    Output: Parser
    '''

    #Define parser
    parser = argparse.ArgumentParser(description="Segmentation Pipeline")

    #Add argumnent to parser
    parser.add_argument("-t","--seg-type", type=str, default = 'single', help="Relative path to folder with stainings")
    parser.add_argument("-fp","--folder-path", type=str, default = 'default', help="Relative path to folder with stainings")
    parser.add_argument("-fm","--file-match", type=str, default = '', help="String to be used to specify files in stainings folder")
    parser.add_argument("-fm2","--file-match2", type=str, default = '', help="String to be used to specify files for second staining in stainins folder")
    parser.add_argument("-nt","--name-tag", type=str, default = 'default', help="Name tag to use for output files. Default is file-match arguments")
    parser.add_argument("-o","--output-folder-path", type=str, default = './processed_data', help="Relative path to folder to save output")
    parser.add_argument("-si","--section-indices", nargs="+", type=str, default=["A1","B1","C1","D1","A2","B2","C2","D2"], help="List of section indices to segment. Default is Resolve naming scheme")
    parser.add_argument("-ic","--index-converter", nargs="+", type=str, default='default', help="List of section indices to use that matches default naming defined in section-index argument")
    parser.add_argument("-ic2","--index-converter2", nargs="+", type=str, default='default', help="List of section indices to use for second staining that matches default naming defined in section-index argument")
    parser.add_argument("-ri","--run-index", type=str, help="List of section indices to segment")
    parser.add_argument("-gf","--grid-fill", type=str, default = 'True', help="Fill grid in stainings before segmentation to avoid slicing of cells. (True/False)")
    parser.add_argument("-d","--diameter", type=float, default = None, help="Diameter to use in cellpose. Default is auto-estimation of diameter")

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
    seg_type = args.seg_type
    path_to_staining_folder = args.folder_path
    file_match = args.file_match
    file_match2 = args.file_match2
    name_tag = args.name_tag
    path_to_output_folder = args.output_folder_path
    section_indices = args.section_indices
    index_converter = args.index_converter
    index_converter2 = args.index_converter2
    run_index = args.run_index
    grid_fill = args.grid_fill
    diameter = args.diameter

    #Define name tag to be used for output files
    if name_tag == 'default':
        name_tag = f"{file_match}_{file_match2}"
        
    #Define path to staining folder using default argument
    if args.folder_path == 'default':
        path_to_staining_folder = f"./raw_data/{run_index}"

    #Dictionary to convert between naming schemes for different section index depending on spatial indices for staining1
    if index_converter == 'default':
        section_index_pairing = dict(zip(section_indices, section_indices))
    else:
        section_index_pairing = dict(zip(section_indices, index_converter))

    #Dictionary to convert between naming schemes for different section index depending on spatial indices for staining2
    if index_converter2 == 'default':
        section_index_pairing2 = dict(zip(section_indices, section_indices))
    else:
        section_index_pairing2 = dict(zip(section_indices, index_converter2))
    
    #Create folder to output ROIs
    os.makedirs(f'{path_to_output_folder}/{run_index}',exist_ok=True)

    #Loop over spatial indices 
    for i, section_index in enumerate(section_indices):

        #Load and prepare image depending of segmentation type
        if seg_type == 'single':

            #Load single staining
            print(f"[{current_time()}] Loading staining using {file_match} as file match for section {section_index} from {path_to_staining_folder} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            image = load_staining(path_to_staining_folder,file_match,section_index,section_index_pairing)

        elif seg_type == "combined":

            #Load staining 1 and staining 2 according to their respective section index
            print(f"[{current_time()}] Loading stainings using {file_match} and {file_match2} as file matches for section {section_index} from {path_to_staining_folder} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            image1 = load_staining(path_to_staining_folder,file_match,section_index,section_index_pairing)
            image2 = load_staining(path_to_staining_folder,file_match2,section_index,section_index_pairing2)

            #Combine stainings by using average pixel value
            print(f"[{current_time()}] Combining stainings for section {section_index} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            image = (image1 + image2)/2

        else:

            print("Segmentation type not supported!")

        if grid_fill == "True":
            #Generate rois for area correction later i spatial pipeline
            print(f"[{current_time()}] Generating rois for area correction {section_index} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            rois_area_correction = generate_rois_for_area_correction(image)

            #Save rois for area correction
            print(f"[{current_time()}] Save rois for area correction {section_index} to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            np.savez(f'{path_to_output_folder}/{run_index}/{section_index}_{name_tag}_area_correction_rois.npz', **rois_area_correction)

            #Fill grids in staining
            print(f"[{current_time()}] Filling grids on section {section_index} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            image = fill_grids(image)

        #Input the combined image to staining_to_rois, which returns cell ROIs for the image
        print(f"[{current_time()}] Running cellpose on section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        rois = staining_to_rois(image,diameter)

        #Save the ROIs
        print(f"[{current_time()}] Saving rois for section {section_index} to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        np.savez(f'{path_to_output_folder}/{run_index}/{section_index}_{name_tag}_cells_rois.npz', **rois)

        #Save the ROIs in ImageJ format
        print(f"[{current_time()}] Saving rois for section {section_index} in ImageJ format to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        save_rois_as_imagej(rois, path_to_output_folder, run_index, section_index, name_tag)


if __name__ == "__main__":
    main()
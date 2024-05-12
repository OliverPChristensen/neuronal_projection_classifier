import os
import sys
from datetime import datetime
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from skimage import io, measure

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
    
    print(f"[{current_time()}] Diameter used for segmentation: {diams}")
    sys.stdout.flush()
    
    #Convert masks from image into ROIs and put them into a dictionary
    rois = {}
    for label in range(1, np.max(masks) + 1):
        mask = (masks == label).astype(np.uint8)
        contour = measure.find_contours(mask, 0.5)
        rois[f"Cell{label}"] = contour[0][:,[1,0]]

    return rois


def save_rois_as_imagej(rois, path_to_output_folder, section_index, name_tag):
    '''
    Input: ROIs in a dictionary, path to folder input, run index, section index

    Converts the ROIs into imageJ format and saves them in a zip file

    Output: NA
    '''

    #Convert ROIs into imageJ format
    rois_imJ = [ImagejRoi.frompoints(outline) for _, outline in rois.items()]

    #Create filename for the ROIs
    file_name = f'{path_to_output_folder}/{section_index}_{name_tag}_imagej_rois.zip'

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
    parser.add_argument("-o","--output-folder-path", type=str, help="Relative path to folder to save output")
    parser.add_argument("-sp","--staining-path", type=str, help="Relative path to staining")
    parser.add_argument("-sp2","--staining-path2", type=str, default = None, help="Relative path to staining 2")
    parser.add_argument("-ri","--run-index", type=str, help="Run index of the section")
    parser.add_argument("-si","--section-index", type=str, help="Section index of the section")
    parser.add_argument("-nt","--name-tag", type=str, default = '', help="Name tag to use for output files")
    parser.add_argument("-gf","--grid-fill", type=str, default = 'True', help="Fill grid in stainings before segmentation to avoid slicing of cells. (True/False)")
    parser.add_argument("-d","--diameter", type=float, default = None, help="Diameter to use in cellpose. Default is auto-estimation of diameter")

    return parser


def main():
    '''
    Input: NA

    Main function that defines parser arguments, loades staining(s), combines them if multiple, performs grid fill if selected and generates image ROIs
    runs cellpose on the (combined) image and creates cell ROIs from the masks outputted, save the cell ROIs in numpy and imageJ format

    Output: NA
    '''

    #Create parser and parse arguments
    parser = create_parser()
    args = parser.parse_args()
    
    #Define variables from parser
    path_to_output_folder = args.output_folder_path
    staining_path = args.staining_path
    staining_path2 = args.staining_path2
    run_index = args.run_index
    section_index = args.section_index
    name_tag = args.name_tag
    grid_fill = args.grid_fill
    diameter = args.diameter
    
    #Create folder to output files
    os.makedirs(f'{path_to_output_folder}',exist_ok=True)

    #Load and prepare image depending of segmentation type
    if staining_path2 == None:

        #Load single staining
        print(f"[{current_time()}] Loading staining using from {staining_path} for section {section_index}")
        sys.stdout.flush()
        image = io.imread(staining_path)

    elif type(staining_path2) == str:

        #Load staining 1 and staining 2 according to their respective section index
        print(f"[{current_time()}] Loading stainings using from {staining_path} and {staining_path2} for section {section_index}")
        sys.stdout.flush()
        image1 = io.imread(staining_path)
        image2 = io.imread(staining_path2)

        #Combine stainings by using average pixel value
        print(f"[{current_time()}] Combining stainings for section {section_index}")
        sys.stdout.flush()
        image = (image1 + image2)/2

    else:

        print("Segmentation type not supported!")

    if grid_fill == "True":
        #Generate rois for area correction later i spatial pipeline
        print(f"[{current_time()}] Generating rois for area correction {section_index}")
        sys.stdout.flush()
        rois_area_correction = generate_rois_for_area_correction(image)

        #Save rois for area correction
        print(f"[{current_time()}] Save rois for area correction {section_index} to {path_to_output_folder}")
        sys.stdout.flush()
        np.savez(f'{path_to_output_folder}/{section_index}_{name_tag}_area_correction_rois.npz', **rois_area_correction)

        #Fill grids in staining
        print(f"[{current_time()}] Filling grids on section {section_index}")
        sys.stdout.flush()
        image = fill_grids(image)

    #Input the combined image to staining_to_rois, which returns cell ROIs for the image
    print(f"[{current_time()}] Running cellpose on section {section_index}")
    sys.stdout.flush()
    rois = staining_to_rois(image,diameter)

    #Save the ROIs
    print(f"[{current_time()}] Saving rois for section {section_index} to {path_to_output_folder}")
    sys.stdout.flush()
    np.savez(f'{path_to_output_folder}/{section_index}_{name_tag}_cells_rois.npz', **rois)

    #Save the ROIs in ImageJ format
    print(f"[{current_time()}] Saving rois for section {section_index} in ImageJ format to {path_to_output_folder}")
    sys.stdout.flush()
    save_rois_as_imagej(rois, path_to_output_folder, section_index, name_tag)


if __name__ == "__main__":
    main()

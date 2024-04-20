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
    file_name = matching_files[np.core.defchararray.find(matching_files,section_index_pairing[section_index]) > 0].item()

    #Define path to staining file
    path_to_file = f"{path_to_staining_folder}/{file_name}"
    
    #Load staining file
    image = io.imread(path_to_file)

    return image


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
    print(f"[{current_time()}] Device is {model.device}")
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
    parser.add_argument("-ic","--index_converter", nargs="+", type=str, default='default', help="List of section indices to use that matches default naming defined in section-index argument")
    parser.add_argument("-ic2","--index_converter2", nargs="+", type=str, default='default', help="List of section indices to use for second staining that matches default naming defined in section-index argument")
    parser.add_argument("-ri","--run-index", type=str, help="List of section indices to segment")
    parser.add_argument("-gf","--grid-fill", type=str, default = 'True', help="Fill grid in stainings before segmentation to avoid slicing of cells. (True/False)")

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

    if name_tag == 'default':
        name_tag = f"{file_match}_{file_match2}"
        
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

            #Load sun1 image
            print(f"[{current_time()}] Loading staining for section {section_index} from {path_to_staining_folder} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            image = load_staining(path_to_staining_folder,file_match,section_index,section_index_pairing)

        elif seg_type == "combined":

            #Load DAPI image and polyT image according to section index
            print(f"[{current_time()}] Loading stainings for section {section_index} from {path_to_staining_folder} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            image1 = load_staining(path_to_staining_folder,file_match,section_index,section_index_pairing)
            image2 = load_staining(path_to_staining_folder,file_match2,section_index,section_index_pairing2)

            #Combine dapi and polyT stainings by using average pixel value
            print(f"[{current_time()}] Combining stainings for section {section_index} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            image = (image1 + image2)/2

        else:

            print("Segmentation type not supported!")

        if grid_fill == "True":
            #Fill grids in staining
            print(f"[{current_time()}] Filling grids on section {section_index} ({i+1}/{len(section_indices)})")
            sys.stdout.flush()
            image = fill_grids(image)

        #Input the combined image to staining_to_rois, which returns cell ROIs for the image
        print(f"[{current_time()}] Running cellpose on section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        rois = staining_to_rois(image)

        #Save the ROIs
        print(f"[{current_time()}] Saving rois for section {section_index} to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        np.savez(f'{path_to_output_folder}/{run_index}/{section_index}_{name_tag}_rois.npz', **rois)

        #Save the ROIs in ImageJ format
        print(f"[{current_time()}] Saving rois for section {section_index} in ImageJ format to {path_to_output_folder}/{run_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        save_rois_as_imagej(rois, path_to_output_folder, run_index, section_index, name_tag)


if __name__ == "__main__":
    main()












def load_dapi_polyT_pair(path_to_staining_folder,dapi_file_match,polyT_file_match,section_index,section_index_pairing):
    '''
    Input: Path to staining folder, DAPI file match string, polyT file match string, section index, dictionary with section index pairing

    Loads DAPI path according to section index, finds corresponding polyT path using section index pairing, loads the DAPI and polyT image

    Output: DAPI image as numpy array, polyT image as numpy array
    '''
    #Load files in folder with stanings
    files = os.listdir(path_to_staining_folder)

    #Fetch all dapi and polyT file names
    dapi_list = np.array([file for file in files if dapi_file_match in file])
    polyT_list = np.array([file for file in files if polyT_file_match in file])

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


"""
#Generate quility report of the segmentation using quality_report
print(f"[{current_time()}] Generating quality report for section {section_index} ({i+1}/{len(section_indices)})")
sys.stdout.flush()
quality_report(image, rois, run_index, section_index)
"""
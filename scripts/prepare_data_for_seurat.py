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
    print(f"{current_time()}: Device is {model.device}")
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

def load_dapi_polyT_pair(path_to_staining_folder,section_index):

    #Dictionary to convert between naming schemes for different sections
    section_index_pairing = {
        "A1": "W0",
        "B1": "W1",
        "C1": "W2",
        "D1": "W3",
        "A2": "W4",
        "B2": "W5",
        "C2": "W6",
        "D2": "W7"
    }

    #Load files in folder with stanings
    files = os.listdir(path_to_staining_folder)

    #Fetch all dapi and polyT file names
    dapi_list = np.array([element for element in files if "DAPI" in element])
    polyT_list = np.array([element for element in files if "Cy5" in element])

    #Find the dapi file name that matches the section index
    dapi_file = dapi_list[np.core.defchararray.find(dapi_list,section_index) > 0].item()

    #Fetch the alternative section index paired to the section index
    alt_section_index = section_index_pairing[section_index]

    #Find the polyT file name that matches the alternative section index
    polyT_file = polyT_list[np.core.defchararray.find(polyT_list, alt_section_index) > 0]

    path_to_dapi_file = f"{path_to_staining_folder}/{dapi_file}"
    path_to_polyT_file = f"{path_to_staining_folder}/{polyT_file}"


    #Load paired dapi and polyT stainings
    dapi_image = io.imread(path_to_dapi_file)
    polyT_image = io.imread(path_to_polyT_file)

    return dapi_file, polyT_file



def main():
    parser = argparse.ArgumentParser(description="Segmentation Pipeline")
    parser.add_argument("-fp","--folder-path", type=str, default = 'default', help="Relative path to folder with stainings")
    parser.add_argument("-o","--output-folder-path", type=str, default = './processed_data', help="Relative path to folder to save output")
    parser.add_argument("-si","--section-indices", nargs="+", type=str, default=["A1","B1","C1","D1","A2","B2","C2","D2"], help="List of section indices to segment")
    parser.add_argument("-ri","--run-index", type=str, help="List of section indices to segment")
    
    args = parser.parse_args()
    
    path_to_output_folder = args.output_folder_path
    section_indices = args.section_indices
    run_index = args.run_index

    if args.folder_path == 'default':
        path_to_staining_folder = f"./raw_data/{run_index}"
    else:
        path_to_staining_folder = args.folder_path


    for i, section_index in enumerate(section_indices):
        dapi_image, polyT_image = load_dapi_polyT_pair(path_to_staining_folder,section_index)

        #Combine dapi and polyT stainings by using average pixel value
        image = (dapi_image + polyT_image)/2

        #Input the combined image to staining_to_rois, which returns cell ROIs for the image
        print(f"{current_time()}: Running cellpose on section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        rois = staining_to_rois(image)

        #Save the ROIs
        print(f"{current_time()}: Saving rois for section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        np.savez(f'{path_to_output_folder}/{run_index}/{section_index}_rois.npz', **rois)

        #Generate quility report of the segmentation using quality_report
        print(f"{current_time()}: Generating quality report for section {section_index} ({i+1}/{len(section_indices)})")
        sys.stdout.flush()
        quality_report(image, rois, run_index, section_index)


if __name__ == "__main__":
    main()



"""


def create_rois(folder_with_stainings,run):
    '''
    Input: Folder with staining images, spatial run index

    Loads polyT and DAPI stainings in corresponding pairs and combines them using the average pixel value
    Forwards the image to stainings_to_rois to generate segmentation
    Saves the ROIs and create quality report of segmentation

    Output: NA
    '''

    #Dictionary to convert between naming schemes for different sections
    name_con_dict = {
        "W0": "A1",
        "W1": "B1",
        "W2": "C1",
        "W3": "D1",
        "W4": "A2",
        "W5": "B2",
        "W6": "C2",
        "W7": "D2"
    }

    #Load files in folder with stanings
    files = os.listdir(folder_with_stainings)

    #Fetch all dapi and polyT file names
    dapi_list = np.array([element for element in files if "DAPI" in element])
    polyT_list = np.array([element for element in files if "Cy5" in element])

    #Define variable for keys in name_con_dict to be used later
    name_con_keys = np.array(list(name_con_dict.keys()))

    #Loop over all polyT file names
    for i, polyT_file in enumerate(polyT_list):

        #Find the section index key present in the polyT file name
        name_con_ele = name_con_keys[np.core.defchararray.find(polyT_file, name_con_keys) > 0]

        #Fetch the value paired to the section index key
        name_con_ele_pair = name_con_dict[name_con_ele.item()]

        #Find the dapi file name that matches the section index value such that the dapi file paired to the polyT file is now found
        dapi_file = dapi_list[np.core.defchararray.find(dapi_list,name_con_ele_pair) > 0][0]

        #Load paired dapi and polyT stainings
        dapi_image = io.imread(f"{folder_with_stainings}/{dapi_file}")
        polyT_image = io.imread(f"{folder_with_stainings}/{polyT_file}")

        #Combine dapi and polyT stainings by using average pixel value
        image = (dapi_image + polyT_image)/2

        #Input the combined image to staining_to_rois, which returns cell ROIs for the image
        print(f"{current_time()}: Running cellpose on section {name_con_ele_pair} ({i+1}/{len(polyT_list)})")
        sys.stdout.flush()
        rois = staining_to_rois(image)

        #Save the ROIs
        print(f"{current_time()}: Saving rois for section {name_con_ele_pair} ({i+1}/{len(polyT_list)})")
        sys.stdout.flush()
        np.savez(f'./processed_data/{run}/{name_con_ele_pair}_rois.npz', **rois)

        #Generate quility report of the segmentation using quality_report
        print(f"{current_time()}: Generating quality report for section {name_con_ele_pair} ({i+1}/{len(polyT_list)})")
        sys.stdout.flush()
        quality_report(image, rois, run, name_con_ele_pair)

    
run = "spatial1"

path_to_staining_folder = "./raw_data/spatial1"

create_rois(folder_with_stainings,run)



"""

test = 2
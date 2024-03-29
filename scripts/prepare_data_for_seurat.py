import os
import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from skimage import io, measure

from shapely.geometry import Polygon

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




def create_rois(folder_with_stainings,run):
    '''
    Input: Folder with staining images, spatial run index

    Loades polyT and dapy stainings in corresponding pairs and combines them using the average pixel value
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
        polyT_image = io.imread(f"{folder_with_stainings}/{dapi_file}")

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

folder_with_stainings = "./raw_data/spatial1"

create_rois(folder_with_stainings,run)

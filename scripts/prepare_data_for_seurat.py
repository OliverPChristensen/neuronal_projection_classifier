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
    cur_time = datetime.now()
    return cur_time.strftime("%Y-%m-%d %H:%M:%S")

def staining_to_rois(image):
    '''
    Input: Folder with polyT and dapi staining

    Loades and combines polyT and dapi staining (average pixel value)
    Runs cellpose on the combined staining image
    Converts masks to list of ROIs

    Output: ROIs for the segemnetation of the stainings
    '''


    model = models.Cellpose(gpu=True, model_type = 'cyto')

    print(f"{current_time()}: Device is {model.device}")
    sys.stdout.flush()
    masks, _, _, diams = model.eval(image, diameter=None, channels=[0,0],
                                         flow_threshold=0.4, do_3D=False)
    
    rois = {}

    for label in range(1, np.max(masks) + 1):
        mask = (masks == label).astype(np.uint8)
        contour = measure.find_contours(mask, 0.5)
        rois[f"Cell{label}"] = contour[0][:,[1,0]]

    return rois

def quality_report(image,rois,image_index,section):

    fig, ax = plt.subplots()

    _ = ax.imshow(image)
    
    for roi in rois:
        #print(type(roi))
        patch = patches.Polygon(rois[roi], closed=True, edgecolor='r', facecolor='none')
        _ = ax.add_patch(patch)
        # Add the rectangle to the plot

    image_split_num = 50

    xdim = int(image.shape[1]/image_split_num)
    ydim = int(image.shape[0]/image_split_num)
    
    os.makedirs(f"./plots/segmentation_report/{image_index}/{section}",exist_ok=True)

    for i in range(image_split_num):
        for j in range(image_split_num):
            
            if np.max(image[ydim*j:ydim*(j + 1),xdim*i:xdim*(i + 1)]) > 0:
                _ = ax.set_xlim(xdim*i, xdim*(i + 1))  # Adjust these values as needed
                _ = ax.set_ylim(ydim*j, ydim*(j + 1))  # Adjust these values as needed

                plt.savefig(f"./plots/segmentation_report/{image_index}/{section}/{section}_quality_report_{i}-{j}.png")




def create_rois(folder_with_stainings,spatial_index):
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

    files = os.listdir(folder_with_stainings)

    dapi_list = np.array([element for element in files if "DAPI" in element])
    polyT_list = np.array([element for element in files if "Cy5" in element])

    name_con_keys = np.array(list(name_con_dict.keys()))

    for i, polyT_file in enumerate(polyT_list):

        name_con_ele = name_con_keys[np.core.defchararray.find(polyT_file, name_con_keys) > 0]

        name_con_ele_pair = name_con_dict[name_con_ele.item()]
        dapi_file = dapi_list[np.core.defchararray.find(dapi_list,name_con_ele_pair) > 0][0]

        dapi_image = io.imread(f"{folder_with_stainings}/{dapi_file}")
        polyT_image = io.imread(f"{folder_with_stainings}/{dapi_file}")

        image = (dapi_image + polyT_image)/2

        print(f"{current_time()}: Running cellpose on image {i+1}/{len(polyT_list)}")
        sys.stdout.flush()
        rois = staining_to_rois(image)

        print(f"{current_time()}: Saving rois for image {i+1}/{len(polyT_list)}")
        sys.stdout.flush()
        np.savez(f'./processed_data/{spatial_index}/{name_con_ele_pair}_rois.npz', **rois)

        print(f"{current_time()}: Generating quality report")
        sys.stdout.flush()
        quality_report(image, rois, spatial_index, name_con_ele_pair)



spatial_index = "spatial1"

folder_with_stainings = "./raw_data/spatial1"

create_rois(folder_with_stainings,spatial_index)

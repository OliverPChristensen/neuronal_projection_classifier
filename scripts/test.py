import skimage as ski
import matplotlib.pyplot as plt
import numpy as np
from skimage import measure
import cv2
from roifile import ImagejRoi, roiwrite
import matplotlib.patches as patches

image = ski.io.imread("./raw_data/spatial1/Panorama_Spatial1_serie2_W0A1_Cy5-class_R11_.tiff")


mask = (image != 0).astype(np.uint8)
contour = measure.find_contours(mask, 0.5)

fig, ax = plt.subplots()

ax.imshow(image)

rois = {}
for i, cont in enumerate(contour):
    rois[str(i)] = cont[:,[1,0]]


np.savez("processed_data/spatial1/test.npz", **rois)


im_crop = image[11000:15000,11000:15000]

plt.imshow(im_crop)
plt.savefig("./plots/pplot.png")


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


im_filled = fill_grids(im_crop)




from cellpose import models

model = models.Cellpose(gpu=True, model_type = 'cyto')


masks, flows, styles, diams = model.eval(im_filled, diameter=None, channels=[0,0],
                                         flow_threshold=0.4, do_3D=False)



from skimage import measure
import matplotlib.patches as patches

#Convert masks from image into ROIs and put them into a dictionary
rois = {}
for label in range(1, np.max(masks) + 1):
    mask = (masks == label).astype(np.uint8)
    contour = measure.find_contours(mask, 0.5)
    rois[f"Cell{label}"] = contour[0][:,[1,0]] + 5000





fig, ax = plt.subplots()

ax.imshow(image)


for _,roi in rois.items():
    patch = patches.Polygon(roi, closed=True, edgecolor='r', facecolor='none')
    ax.add_patch(patch)


plt.savefig("./plots/pplot.png")





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

import os

save_rois_as_imagej(rois, "./processed_data","spatial1_gridfilled", "test2")






a = np.array([1,1])
os.makedirs("./processed_data/test",exist_ok=True)
np.savez('./processed_data/test/test.npz', a)


image = ski.io.imread('./raw_data/spatial1/Panorama_Spatial1_serie2_W3A1_Cy5-class_R11_.tiff')
#image = io.imread('./raw_data/Panorama_EXP-2_W1A1_Cy5-class_R9_.tiff')
#ski.io.imshow(crop_im)

crop_im = image[0:2000,0:2000]

blurred_image = ski.filters.gaussian(crop_im, sigma=8) 

fig, ax = plt.subplots()

xdim = int(image.shape[0])
ydim = int(image.shape[1])


_ = ax.set_xlim(0, 1000)  # Adjust these values as needed

ax.imshow(image)
plt.savefig("./plots/test2.png")

np.max(images[0].get_array())
images = ax.get_images()
ax.axis('off')  # Turn off axis
plt.show()


plt.savefig("./plots/test2.png")

crop_im.sum()
ski.io.imsave("./plots/test2.png",image)

from cellpose import models

model = models.Cellpose(gpu=True, model_type = 'cyto')


masks, flows, styles, diams = model.eval(crop_im, diameter=None, channels=[0,0],
                                         flow_threshold=0.4, do_3D=False)




np.unique(masks)




import numpy as np
import skimage
from skimage import measure
import matplotlib.patches as patches

# Find contours for each label separately
contour_image = np.zeros_like(masks)

rois = []

for label in range(1, np.max(masks) + 1):
    mask = (masks == label).astype(np.uint8)
    contour = measure.find_contours(mask, 0.5)
    rois.append(contour[0][:,[1,0]])



crop_im.shape[0]/5

type(rois[0][0])
contour_image.shape

rois[0][0][:,[0,1]]

roi_x = rois[0][0][:,0]
roi_y = rois[0][0][:,1]
roi_width = np.max(roi_x) - np.min(roi_x)
roi_height = np.max(roi_y) - np.min(roi_y)

for roi in rois:
    print(type(roi))
    patch = patches.Polygon(roi, closed=True, edgecolor='r', facecolor='none')
    ax.add_patch(patch)
# Add the rectangle to the plot


ax.set_xlim(0, 1000)  # Adjust these values as needed
ax.set_ylim(0, 1000)  # Adjust these values as needed

plt.savefig("./plots/test2.png")








model.device

parameter_device = next(model.parameters()).device


import torch


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print("Model is on:", parameter_device)

a = 2

import matplotlib.pyplot as plt
import numpy as np

# Generate data
x = np.linspace(0, 10, 100)  # 100 points between 0 and 10
y = np.sin(x)  # Compute sine values for each x

# Plot
plt.plot(x, y)
plt.title('Sine Wave')
plt.xlabel('x')
plt.ylabel('sin(x)')
plt.grid(True)
plt.show()

from shapely.geometry import Polygon, LineString, MultiPolygon

# Convert NumPy arrays to Shapely polygons
polygons = []
for roi in rois:
    polygon = Polygon(roi)
    polygons.append(polygon)


# Define the main geometric object (e.g., polygon, line, etc.)
main_object = Polygon([(0, 0), (0, 2), (2, 2), (2, 0)])

# Define a list of other geometric objects
other_objects = [
    Polygon([(1, 1), (1, 3), (3, 3), (3, 1)]),
    Polygon([(2, -1), (2, 1), (4, 1), (4, -1)]),
    LineString([(0, 1), (2, 1)])
]

# TODO: Putte objecter i dinctionaries, så deres key er navn på enten celle eller typen af transcript?
# burde være fint fordi transcirpter kan bare tilføjes til seurat og deres association med celle sker via count matrix. Og seg og centroids er kun basseret på rois

# Find intersected objects
intersected_objects = []
for obj in other_objects:
    if main_object.intersects(obj):
        intersected_objects.append(obj)

# If you want to visualize the intersected objects
# You can plot them using Matplotlib or any other plotting library

# Print intersected objects (optional)
print("Intersected objects:")
for obj in intersected_objects:
    print(obj)




from shapely.geometry import Polygon, MultiPolygon

# Create some polygons
polygon1 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
polygon2 = Polygon([(1.5, 0), (2.5, 0), (2.5, 1), (1.5, 1)])

# Create a Multipolygon
multipolygon = MultiPolygon([polygon1, polygon2])

# Create another polygon
polygon3 = Polygon([(0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (0.25, 0.75)])

# Check intersection between the Multipolygon and polygon3
intersection = multipolygon.intersection(polygon3)

print(intersection)






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

## Unused functions
from shapely.geometry import Polygon, Point
import pandas as pd


def transcripts_geometry(trascript_path):

    transcripts_df = pd.read_csv('./raw_data/EXP-2_B1-1_results.txt', sep='\t', header=None)
    transcripts_df.columns = ['x','y','UMI','ID','NA']

    transcripts_IDs = transcripts_df['ID'].unique()
    transcripts_IDs = transcripts_IDs[['FP ' not in item for item in transcripts_IDs]]

    # Create a list to store Shapely Point objects
    transcripts_list = {}

    # Iterate over the DataFrame and create Point objects
    for i, transcripts_ID in enumerate(transcripts_IDs):

        point_list = []
        transcripts_subset = transcripts_df[transcripts_df['ID'] == transcripts_ID]

        for index, row in transcripts_subset.iterrows():
            point = Point(row['x'], row['y'])
            point_list.append(point)

        
        print(f"trancsript {i}/100")
        

        transcripts_list[transcripts_ID] = point_list

        
    
    
    return transcripts

def create_count_matrix(rois, transcripts_path):

    transcripts_df = pd.read_csv('./raw_data/EXP-2_B1-1_results.txt', sep='\t', header=None)
    transcripts_df.columns = ['x','y','UMI','ID','NA']

    transcripts_IDs = transcripts_df['ID'].unique()
    transcripts_IDs = transcripts_IDs[['FP ' not in item for item in transcripts_IDs]]
    
    count_matrix = pd.DataFrame(0, index=rois.keys(), columns=transcripts_IDs)
    for roi in rois:
        polygon = Polygon(rois[roi])
        # Iterate over the DataFrame and create Point objects
        for i, transcripts_ID in enumerate(transcripts_IDs):

            transcripts_subset = transcripts_df[transcripts_df['ID'] == transcripts_ID]

            count = 0

            for index, row in transcripts_subset.iterrows():
                point = Point(row['x'], row['y'])

                if polygon.intersection(point):
                    count += 1
                
            
            count_matrix.loc[roi,transcripts_ID]

            
            print(f"trancsript {i}/100")


    return
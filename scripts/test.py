import skimage as ski
import matplotlib.pyplot as plt
from cellpose import io
import numpy as np

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

model = models.Cellpose(gpu=False, model_type = 'cyto')


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
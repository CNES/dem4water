import matplotlib.pyplot as plt

# from skimage.segmentation import morphological_chan_vese
import numpy as np
import rasterio as rio
from rasterio.plot import reshape_as_image, reshape_as_raster
from skimage import color, io, measure

with rio.open(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Naussac/dem_extract_Naussac.tif"
) as im:
    image = im.read()
    # image = morphological_chan_vese(image, num_iter=5, init_level_set="checkerboard")
    # with rio.open(
    #     "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/snake.tif",
    #     "w",
    #     im.profile,
    # ) as out:
    #     out.write(image)


# Construct some test data
# x, y = np.ogrid[-np.pi:np.pi:100j, -np.pi:np.pi:100j]
# r = np.sin(np.exp(np.sin(x)**3 + np.cos(y)**2))
print(image.shape)
r = reshape_as_image(image)
r = r.reshape(r.shape[0], r.shape[1])

print(r.shape)
# Find contours at a constant value of 0.8
contours = measure.find_contours(r, 943)
print(contours)
# Display the image and plot all contours found
fig, ax = plt.subplots()
ax.imshow(r, cmap=plt.cm.gray)

for contour in contours:
    ax.plot(contour[:, 1], contour[:, 0], "r", linewidth=2)

ax.axis("image")
ax.set_xticks([])
ax.set_yticks([])
plt.show()

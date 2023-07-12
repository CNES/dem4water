import numpy as np
import rasterio as rio
import scipy


def compute_gradient_product(water_binary_mask, dem_extract, output_raster):
    """."""

    with rio.open(water_binary_mask, "r") as wmap:
        with rio.open(dem_extract, "r") as dem:
            wmap_array = wmap.read(1)
            dem_array = dem.read(1)

            mask_lisse = scipy.ndimage.gaussian_filter(wmap_array, sigma=5)
            dx_mask_sobel = scipy.ndimage.sobel(mask_lisse, axis=0, mode="constant")
            dy_mask_sobel = scipy.ndimage.sobel(mask_lisse, axis=1, mode="constant")

            image_lisse = scipy.ndimage.gaussian_filter(dem_array, sigma=5)

            # apply Sobel operator to get the gradient of the DEM image
            dx_sobel = scipy.ndimage.sobel(image_lisse, axis=0, mode="constant")
            dy_sobel = scipy.ndimage.sobel(image_lisse, axis=1, mode="constant")
            # compute pixel-wise dot product between the two vector fields
            # (normal contour vectors and gradient of the DEM)
            im_produit_sobel = dx_sobel * dx_mask_sobel + dy_sobel * dy_mask_sobel

            # set interior of the contour to 0
            im_produit_sobel = np.where(wmap_array == 0, im_produit_sobel, 0)
            # binarise to keep only descending slope (dam output)
            im_produit_sobel = np.where(im_produit_sobel > 0, 1, 0)
            profile = wmap.profile
            profile.update(dtype=rio.float64, count=1, compress="lzw")
            with rio.open(output_raster, "w", **profile) as output:
                output.write(im_produit_sobel.astype(rio.float64), 1)


compute_gradient_product(
    "/home/btardy/Documents/activites/WATER/extract_montbel/wmap_extract-Montbel.tif",
    "/home/btardy/Documents/activites/WATER/extract_montbel/dem_extract-Montbel.tif",
    "/home/btardy/Documents/activites/WATER/extract_montbel/test_gradient_product_bin.tif",
)

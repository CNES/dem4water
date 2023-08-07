# import cv2
import numpy as np
import rasterio as rio
import scipy

# def get_rotated(array):
#     """
#     Function computing the principal axis of a binary mask and turning the mask to align the principal axis with the y-axis. The principal axis corresponds to the direction of maximum variance
#     i.e. the direction associated with the highest eigenvalue of the covariance matrix.

#     Arguments:

#         array (array): binary mask

#     Returns:

#         rotated_array (array): rotated binary mask

#         theta (float): angle (in rad) between the principal axis and the y-axis of the image

#     """

#     y, x = np.nonzero(array)
#     x = x - np.mean(x)
#     y = y - np.mean(y)
#     coords = np.vstack([x, y])
#     cov = np.cov(coords)
#     evals, evecs = np.linalg.eig(cov)
#     sort_indices = np.argsort(evals)[::-1]
#     x_v1, y_v1 = evecs[:, sort_indices[0]]  # Eigenvector with largest eigenvalue
#     x_v2, y_v2 = evecs[:, sort_indices[1]]
#     theta = np.arctan((x_v1) / (y_v1))
#     rotation_mat = np.matrix(
#         [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
#     )
#     transformed_mat = rotation_mat * coords
#     rotated_array = scipy.ndimage.rotate(array, 360 - np.rad2deg(theta))

#     return rotated_array, theta


def compute_gradient_product(water_binary_mask, dem_extract, output_raster):
    """."""

    with rio.open(water_binary_mask, "r", dtype=rio.float64) as wmap:
        with rio.open(dem_extract, "r", dtype=rio.float64) as dem:
            wmap_array = wmap.read(1)
            dem_array = dem.read(1)
            # print(wmap_array.dtype, dem_array.dtype)
            # dem_array = cv2.resize(
            #     dem_array.astype(float), (250, 250), interpolation=cv2.INTER_LINEAR
            # )
            wmap_array = np.where(wmap_array > 0.1, 255, 0)
            # wmap_array = cv2.resize(
            #     wmap_array.astype(float), (250, 250), interpolation=cv2.INTER_LINEAR
            # )
            profile = wmap.profile
            with rio.open(
                output_raster.replace(".tif", "bin_wmap.tif"), "w", **profile
            ) as output:
                output.write(wmap_array.astype(rio.float64), 1)
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
            # im_produit_sobel = scipy.ndimage.rotate(
            #     im_produit_sobel,
            #     360 - np.rad2deg(get_rotated(wmap_array)[1]),
            #     reshape=False,
            #     order=0,
            # )
            profile.update(dtype=rio.float64, count=1, compress="lzw")
            with rio.open(output_raster, "w", **profile) as output:
                output.write(im_produit_sobel.astype(rio.float64), 1)


compute_gradient_product(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/extract_marne_surfwater_no_holes.tif",
    "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/dem_extract_Marne-Giffaumont.tif",
    "/home/btardy/Documents/activites/WATER/GDP/tests_bds/marne_giffaumont_bin_surfwater_no_holes.tif",
)
# compute_gradient_product(
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/extract_marne_grand.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/dem_extract_Marne-Giffaumont.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/tests_bds/marne_giffaumont_bin_grand.tif",
# )
# compute_gradient_product(
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/extract_marne_pekel.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/dem_extract_Marne-Giffaumont.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/tests_bds/marne_giffaumont_bin_pekel.tif",
# )
# compute_gradient_product(
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/extract_marne_surf.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/dem_extract_Marne-Giffaumont.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/tests_bds/marne_giffaumont_bin_surf.tif",
# )

# compute_gradient_product(
#     "/home/btardy/Documents/activites/WATER/GDP/bin_wholes/Laparan.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Laparan/dem_extract_Laparan.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/tests_wholes/laparan_bin.tif",
# )
# compute_gradient_product(
#     "/home/btardy/Documents/activites/WATER/GDP/bin_wholes/Marne.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/dem_extract_Marne-Giffaumont.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/tests_wholes/marne_giffaumont_bin.tif",
# )
# compute_gradient_product(
#     "/home/btardy/Documents/activites/WATER/GDP/bin_wholes/Naussac.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Naussac/dem_extract_Naussac.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/tests_wholes/naussac_bin.tif",
# )
# compute_gradient_product(
#     "/home/btardy/Documents/activites/WATER/GDP/bin_wholes/Vaufrey.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Vaufrey/dem_extract_Vaufrey.tif",
#     "/home/btardy/Documents/activites/WATER/GDP/tests_wholes/vaufrey_bin.tif",
# )

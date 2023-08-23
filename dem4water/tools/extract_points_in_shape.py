import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio import Affine
from rasterio.mask import mask
from shapely.geometry import Point, mapping

shapefile = gpd.read_file(
    "/home/btardy/Documents/activites/WATER/GDP/test_chain_3/cutline.geojson"
)
shapefile.geometry = shapefile.geometry.buffer(300)
geoms = shapefile.geometry.values

with rasterio.open(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/dem_extract_Marne-Giffaumont.tif"
) as src:
    out_image, out_transform = mask(src, geoms, crop=True)

# no data values of the original raster
no_data = src.nodata
print(no_data)

# extract the values of the masked array
data = out_image[0, :, :]
# extract the row, columns of the valid values
row, col = np.where(data != no_data)
rou = np.extract(data != no_data, data)

# affine import Affine
T1 = out_transform * Affine.translation(0.5, 0.5)  # reference the pixel centre
rc2xy = lambda r, c: (c, r) * T1

d = pd.DataFrame({"col": col, "row": row, "ROU": rou})
# coordinate transformation
d["x"] = d.apply(lambda row: rc2xy(row.row, row.col)[0], axis=1)
d["y"] = d.apply(lambda row: rc2xy(row.row, row.col)[1], axis=1)
# geometry

gdf = gpd.GeoDataFrame(
    d,
    geometry=d.apply(lambda row: Point(row["x"], row["y"]), axis=1),
    crs=shapefile.crs,
)
print(gdf)
# save to a shapefile
gdf.to_file("/home/btardy/Documents/activites/WATER/GDP/test_chain_3/mnt.geojson")

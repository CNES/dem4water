import geopandas as gpd
import pandas as pd


def poly_to_points(feature, poly):
    return {feature: poly.exterior.coords}


algo = "inpe"
gdp_vector_file = f"/home/btardy/Documents/activites/WATER/GDP/tests_bds/marne_giffaumont_bin_{algo}_no_holes.geojson"
water_body_vector = (
    "/home/btardy/Documents/activites/WATER/GDP/extract"
    f"/Marne-Giffaumont/extract_marne_{algo}_noholes.geojson"
)

gdf_wb = gpd.GeoDataFrame().from_file(water_body_vector)

gdf_wb = gdf_wb.to_crs(2154)
gdf_wb = gdf_wb.explode(ignore_index=True)
gdf_wb["points"] = gdf_wb.apply(
    lambda p: poly_to_points(p["id"], p["geometry"]), axis=1
)

points = list(gdf_wb.points)

ident = []
coordx = []
coordy = []
for i, poly in enumerate(points):
    for _, p in points[i].items():
        for p1 in p:
            print(p1[0], p1[1])
            ident.append(i)
            coordx.append(p1[0])
            coordy.append(p1[1])
df = pd.DataFrame()
df["x"] = coordx
df["y"] = coordy
df["ident"] = ident
df["fidx"] = range(len(df.index))

gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x, df.y), crs=gdf_wb.crs)
gdf.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/test_python/test_contour_wb.geojson"
)

import geopandas as gpd
import pandas as pd


def linestring_to_points(feature, line):
    return {feature: line.coords}


def poly_to_points(feature, poly):
    return {feature: poly.exterior.coords}


algo = "inpe"
gdp_vector_file = f"/home/btardy/Documents/activites/WATER/GDP/tests_bds/marne_giffaumont_bin_{algo}_no_holes.geojson"

gdf_gdp = gpd.GeoDataFrame().from_file(gdp_vector_file)

gdf_gdp = gdf_gdp.loc[gdf_gdp.DN == 1]
gdf_gdp = gdf_gdp.to_crs(2154)
gdf_gdp = gdf_gdp.loc[gdf_gdp.geometry.area > 2000]

gdf_gdp["geometry"] = gdf_gdp.geometry.buffer(-20)
gdf_gdp = gdf_gdp.explode(ignore_index=True)
gdf_gdp["geometry"] = gdf_gdp.geometry.simplify(15)
# lines_gdp = gpd.GeoDataFrame(
#     {"ident": range(len(gdf_gdp.index))}, geometry=gdf_gdp.boundary, crs=2154
# )

# lines_gdp["points"] = lines_gdp.apply(
#     lambda x: [y for y in x["geometry"].coords], axis=1
# )
gdf_gdp["points"] = gdf_gdp.apply(
    lambda p: poly_to_points(p["DN"], p["geometry"]), axis=1
)
points = list(gdf_gdp.points)

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

gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x, df.y), crs=gdf_gdp.crs)
gdf.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/test_python/test_poly_to_points.geojson"
)

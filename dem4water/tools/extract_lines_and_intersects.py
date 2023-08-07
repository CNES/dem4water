import geopandas as gpd

algo = "inpe"
gdp_vector_file = f"/home/btardy/Documents/activites/WATER/GDP/tests_bds/marne_giffaumont_bin_{algo}_no_holes.geojson"
water_body_vector = (
    "/home/btardy/Documents/activites/WATER/GDP/extract"
    f"/Marne-Giffaumont/extract_marne_{algo}_noholes.geojson"
)

gdf_gdp = gpd.GeoDataFrame().from_file(gdp_vector_file)
gdf_wb = gpd.GeoDataFrame().from_file(water_body_vector)
gdf_gdp = gdf_gdp.to_crs(2154)
gdf_gdp = gdf_gdp.loc[gdf_gdp.DN == 1]
gdf_wb = gdf_wb.to_crs(2154)

gdf_ori = gpd.GeoDataFrame().from_file(water_body_vector)
gdf_ori = gdf_wb.to_crs(2154)
gdf_gdp["geometry"] = gdf_gdp.geometry.buffer(5)
gdf_wb["geometry"] = gdf_wb.geometry.buffer(15)

intersect = gdf_wb.overlay(gdf_gdp, how="intersection")
cols = [x for x in intersect.columns if not "geometry" in x]

intersect["original_label"] = range(len(intersect.index))
intersect = intersect.drop(columns=cols)

intersect.to_file(
    f"/home/btardy/Documents/activites/WATER/GDP/test_python/{algo}/test_intersect.geojson"
)

polys = intersect.explode(ignore_index=True)
# minimal object size set to 10 pixels
# polys = polys.loc[polys.geometry.area > 1000]
polys["indent"] = range(len(polys.index))
polys.to_file(
    f"/home/btardy/Documents/activites/WATER/GDP/test_python/{algo}/test_explode.geojson"
)
print(polys.boundary)
# polys to line
lines_wb = gpd.GeoDataFrame(
    {"idx": range(len(gdf_ori.index))}, geometry=gdf_ori.boundary, crs=2154
)
lines_gdp = gpd.GeoDataFrame(
    {"ident": range(len(polys.index))}, geometry=polys.boundary, crs=2154
)

lines_wb.to_file(
    f"/home/btardy/Documents/activites/WATER/GDP/test_python/{algo}/test_lines_wb.geojson"
)
lines_gdp.to_file(
    f"/home/btardy/Documents/activites/WATER/GDP/test_python/{algo}/test_lines_gdp.geojson"
)

# intersect lines ?
intersect_lines = lines_wb.overlay(lines_gdp, how="intersection", keep_geom_type=False)
intersect_lines = intersect_lines.explode(ignore_index=True)
intersect_lines.to_file(
    f"/home/btardy/Documents/activites/WATER/GDP/test_python/{algo}/test_intersect_lines.geojson"
)

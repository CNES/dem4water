import geopandas as gpd
from shapely.geometry import LineString

bd_file = "/home/btardy/Documents/activites/WATER/Barrages-15m-France_INPE-V0_340barrages_20221002/392-retenues-pourLoiZSV_V4_sans_tampon.geojson"
gdf = gpd.read_file(bd_file)
gdf_w = gdf[gdf.DAM_NAME == "Alesani"]
gdf_w = gdf_w.to_crs(2154)

list_of_points = list(gdf_w.geometry.values[0].exterior.coords)
distance = 10
new_list_of_points = []
while len(list_of_points) >= 2:
    ori_point = list_of_points.pop(0)  # [index]
    new_list_of_points.append(ori_point)
    # if index + 1 == len(list_of_points):
    #     print("loop ended")
    #     break

    target_point = list_of_points.pop(0)
    line = LineString([ori_point, target_point])
    if distance < line.length:
        print("Too spaced generate point", line.length)
        new_point = line.interpolate(distance)
        list_of_points = [new_point, target_point] + list_of_points
        # new_list_of_points.append(new_point)
    else:
        print("Close enough", line.length)
        # new_list_of_points.append(target_point)
        list_of_points = [target_point] + list_of_points
print(new_list_of_points[-1])
print(new_list_of_points[0])
print(list_of_points)
input("wait")
gdf_w.geometry = [LineString(new_list_of_points)]
gdf_w.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/test_sampling_points_alesani.geojson"
)

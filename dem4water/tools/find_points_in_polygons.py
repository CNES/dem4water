import math

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point


def autoroute_points_df(points_df, x_col="e", y_col="n"):
    """
    Function, that converts a list of random points into ordered points, searching for the shortest possible distance between the points.
    Author: Marjan Moderc, 2016
    """
    points_list = points_df[[x_col, y_col]].values.tolist()

    # arrange points in by ascending Y or X
    points_we = sorted(points_list, key=lambda x: x[0])
    points_sn = sorted(points_list, key=lambda x: x[1])

    # Calculate the general direction of points (North-South or West-East) - In order to decide where to start the path!
    westmost_point = points_we[0]
    eastmost_point = points_we[-1]

    deltay = eastmost_point[1] - westmost_point[1]
    deltax = eastmost_point[0] - westmost_point[0]
    alfa = math.degrees(math.atan2(deltay, deltax))
    azimut = (90 - alfa) % 360

    # If main directon is towards east (45°-135°), take westmost point as starting line.
    if azimut > 45 and azimut < 135:
        points_list = points_we
    elif azimut > 180:
        raise Exception(
            "Error while computing the azimuth! It cant be bigger then 180 since first point is west and second is east."
        )
    else:
        points_list = points_sn

    # Create output (ordered df) and populate it with the first one already.
    ordered_points_df = pd.DataFrame(columns=points_df.columns)
    ordered_points_df = ordered_points_df.append(
        points_df.loc[
            (points_df[x_col] == points_list[0][0])
            & (points_df[y_col] == points_list[0][1])
        ]
    )

    for iteration in range(0, len(points_list) - 1):
        already_ordered = ordered_points_df[[x_col, y_col]].values.tolist()

        current_point = already_ordered[-1]  # current point
        possible_candidates = [
            i for i in points_list if i not in already_ordered
        ]  # list of candidates

        distance = 10000000000000000000000
        best_candidate = None
        for candidate in possible_candidates:
            current_distance = Point(current_point).distance(Point(candidate))
            if current_distance < distance:
                best_candidate = candidate
                distance = current_distance

        ordered_points_df = ordered_points_df.append(
            points_df.loc[
                (points_df[x_col] == best_candidate[0])
                & (points_df[y_col] == best_candidate[1])
            ]
        )

    return ordered_points_df


in_file = (
    "/home/btardy/Documents/activites/WATER/GDP/tests_bds/inpe_no_holes_buff.geojson"
)

gdf = gpd.GeoDataFrame.from_file(in_file)
gdf = gdf.explode(index_parts=False)
gdf["area"] = gdf.geometry.area
max_area = np.max(list(gdf.area))
print(max_area)
gdf = gdf.loc[gdf.area > 0.05 * max_area]
gdf.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/tests_bds/inpe_wsmallarea.geojson"
)
print(gdf)
ser = gpd.GeoSeries(gdf.geometry)
gdf2 = gpd.GeoDataFrame(geometry=ser.representative_point())
gdf2["ident"] = np.ones(len(gdf2.index))
df = pd.DataFrame({"e": gdf2.geometry.x, "n": gdf2.geometry.y})

print(df)
print(gdf2)
df = autoroute_points_df(df, x_col="e", y_col="n")
gdf2.geometry.x = df.e
gdf2.geometry.y = df.n

gdf2.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/tests_bds/rep_points_sort.geojson"
)
print(gdf2)
gdf2.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/tests_bds/rep_points_marne_tampon15m_no_holes.geojson"
)
# convert to line
# treat each `ID` group of points as a line
lines = gdf2.groupby(["ident"])["geometry"].apply(lambda x: LineString(x.tolist()))

# store as a GeodataFrame and add 'ID' as a column (currently stored as the 'index')
lines = gpd.GeoDataFrame(lines, geometry="geometry", crs=gdf.crs)
lines.reset_index(inplace=True)
lines.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/tests_bds/extract_marne_tampon15m_cutline_no_holes.geojson"
)

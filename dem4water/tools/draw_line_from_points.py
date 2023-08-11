import geopandas as gpd
import numpy as np
import pandas as pd

# from scipy.sparse import csr_matrix
# from scipy.sparse.csgraph import shortest_path


# print(gdf)


def run_left(dist_mat, points_to_view, ordored_points, start_point, max_val):
    """."""
    run_stop = False
    current_point = start_point
    while points_to_view and not run_stop:
        points_to_view.remove(current_point)
        ind_min = np.argmin(dist_mat[current_point], axis=0)
        min_val = dist_mat[current_point][ind_min]
        if ind_min == current_point + 1:
            print("je suis continu")
            # j'interdis que cette paire soit à nouveau visible
            dist_mat[current_point][:] = max_val + 10
            dist_mat[:, current_point] = max_val + 10
            dist_mat[:, ind_min] = max_val + 10
            ordored_points = [ind_min] + ordored_points
            current_point = ind_min
        else:
            run_stop = True
    return dist_mat, points_to_view, ordored_points


def run_right():
    """."""


def run_cutline():
    gdf = gpd.GeoDataFrame().from_file(
        "/home/btardy/Documents/activites/WATER/GDP/test_chain/filterd_points_no_dup.geojson"
    )
    for ident in [17]:  # gdf.ident.unique():
        gdf_work = gdf[gdf.ident == ident]
        gdf_work = gdf_work.reset_index(drop=True)
        gdf_res = gdf_work.geometry.apply(lambda g: gdf_work.distance(g))
        max_val = max(gdf_res.max())
        max_allowed_dist = max_val / 10.0
        print(max_allowed_dist)
        for i in range(len(gdf_res.index)):
            gdf_res[i][i] = max_val + 1
        array = gdf_res.values.copy()
        ori_array = gdf_res.values.copy()
        points_to_view = list(gdf_res.index)
        ordored_points = [points_to_view[0]]
        current_point = points_to_view[0]
        # # init
        # array, points_to_view, ordored_points = run_left(
        #     array, points_to_view, ordored_points, current_point, max_val
        # )

        previous_right = False
        previous_left = False
        while points_to_view:
            # array, points_to_view, ordored_points = run_left(
            #     array, points_to_view, ordored_points, current_point, max_val
            # )
            # print(points_to_view)
            # input("wait")
            # print(array)
            # print(ori_array)
            if current_point in points_to_view:
                points_to_view.remove(current_point)
            print(points_to_view)
            print("ord : ", ordored_points)
            print("curr point: ", current_point)
            # # je cherche le min du point courant
            ind_min = np.argmin(array[current_point], axis=0)
            # # print(array[current_point])
            # print("curr point : ", current_point)
            print("indice min : ", ind_min)
            min_val = array[current_point][ind_min]
            print("distance : ", min_val)
            if ind_min == current_point + 1:
                print("je suis continu à gauche")
                previous_left = True
                previous_right = False
                # j'interdis que cette paire soit à nouveau visible
                array[current_point][:] = max_val + 10
                array[:, current_point] = max_val + 10
                array[:, ind_min] = max_val + 10
                ordored_points = [ind_min] + ordored_points
                current_point = ind_min
            elif ind_min == current_point - 1:
                print("je suis continu à droite")
                previous_left = False
                previous_right = True
                # j'interdis que cette paire soit à nouveau visible
                array[current_point][:] = max_val + 10
                array[:, current_point] = max_val + 10
                array[:, ind_min] = max_val + 10
                ordored_points.append(ind_min)
                current_point = ind_min

            else:
                print("je ne suis pas continu")
                candidat = ordored_points[-1]
                right_dist = ori_array[candidat][ind_min]
                print("candidat : ", candidat)
                print("max dist : ", max_allowed_dist)
                print("right_dist: ", right_dist)
                print("min_val", min_val)
                if right_dist < min_val and right_dist < max_allowed_dist:
                    print("le point à droite est plus proche")
                    previous_left = False
                    array[candidat][:] = max_val
                    array[:, candidat] = max_val
                    array[:, ind_min] = max_val
                    ordored_points.append(ind_min)
                    # points_to_view.remove(ind_min)
                    current_point = ind_min
                else:
                    print("change de sens")
                    if previous_right:
                        # go to left
                        current_point = ordored_points[0]

                    elif previous_left:
                        # go to right
                        current_point = ordored_points[-1]
                    else:
                        print("Coupe la ligne")

            # input("wait")

    # graph = csr_matrix(array)
    # print(graph)
    # dist_matrix, predecessors = shortest_path(
    #     csgraph=graph, directed=True, indices=0, return_predecessors=True
    # )
    # print(dist_matrix.shape)

    # print("=" * 10)
    # print("\n" * 4)
    # print(predecessors)
    # ind = np.unravel_index(np.argmax(array, axis=None), array.shape)

    # print(ind)
    # print(gdf_res.idxmax()[48])
    # print(gdf_res.idxmax().idxmax())
    # print(max_val)
    # print(gdf_res)
    # df_temp = gdf_res.sort_values([0], ascending=True)
    # print(df_temp)
    # print(gdf_res[0][17])
    # print(gdf_res[61][62])
    #     view_id = [0]


#     gdf1 = gdf_work.iloc[[0]]
#     gdf2 = gdf_work.iloc[1:]
#     print(gdf1)
#     print(gdf2)
#     dist = gpd.sjoin_nearest(gdf_work, gdf_work, distance_col="distance")
#     # new_id = dist["idx_right"].values[0]
#     # print(new_id)
#     print(dist)

# random coordinates
# gdf_1 = gpd.GeoDataFrame(geometry=gpd.points_from_xy([0, 0, 0], [0, 90, 120]))
# gdf_2 = gpd.GeoDataFrame(geometry=gpd.points_from_xy([0, 0], [0, -90]))
run_cutline()

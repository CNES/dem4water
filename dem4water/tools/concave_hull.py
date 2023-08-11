import geopandas as gpd
import numpy as np
import pygeos
from scipy.spatial import Delaunay
from shapely.geometry import MultiLineString
from shapely.ops import cascaded_union, polygonize


def concave_hull(points_gdf, alpha=35):
    """
    Compute the concave hull (alpha shape) of a GeoDataFrame of points.


    Parameters
    ==========
    points_gdf : gpd.GeoDataFrame
      GeoDataFrame of points.

    alpha: int

      alpha value to influence the gooeyness of the border. Smaller numbers
      don't fall inward as much as larger numbers. Too large, and you lose everything!
    """
    if len(points_gdf) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return points_gdf.unary_union.convex_hull

    coords = pygeos.coordinates.get_coordinates(points_gdf.geometry.values.data)
    tri = Delaunay(coords)
    triangles = coords[tri.vertices]
    a = (
        (triangles[:, 0, 0] - triangles[:, 1, 0]) ** 2
        + (triangles[:, 0, 1] - triangles[:, 1, 1]) ** 2
    ) ** 0.5
    b = (
        (triangles[:, 1, 0] - triangles[:, 2, 0]) ** 2
        + (triangles[:, 1, 1] - triangles[:, 2, 1]) ** 2
    ) ** 0.5
    c = (
        (triangles[:, 2, 0] - triangles[:, 0, 0]) ** 2
        + (triangles[:, 2, 1] - triangles[:, 0, 1]) ** 2
    ) ** 0.5
    s = (a + b + c) / 2.0
    areas = (s * (s - a) * (s - b) * (s - c)) ** 0.5
    circums = a * b * c / (4.0 * areas)
    filtered = triangles[circums < (1.0 / alpha)]
    edge1 = filtered[:, (0, 1)]
    edge2 = filtered[:, (1, 2)]
    edge3 = filtered[:, (2, 0)]
    edge_points = np.unique(np.concatenate((edge1, edge2, edge3)), axis=0).tolist()
    m = MultiLineString(edge_points)
    print(m)
    triangles = list(polygonize(m))
    return gpd.GeoDataFrame(
        {"geometry": [cascaded_union(triangles)]}, index=[0], crs=points_gdf.crs
    )

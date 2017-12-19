"""

"""

import pandas as pn
from .kml_to_plane import *
from .struct_geo import *

# TODO: we should rework the storage of xyz coordinates from single .x, .y, .z class-variables to np.ndarray


def dip(normal_vec):
    return np.arccos(normal_vec[2]) / np.pi * 180.


def dip_dir(normal_vec):
    # +/+
    if normal_vec[0] >= 0 and normal_vec[1] > 0:
        return np.arctan(normal_vec[0]/normal_vec[1]) / np.pi * 180.
    # border cases where arctan not defined:
    elif normal_vec[0] > 0 and normal_vec[1] == 0:
        return 90
    elif normal_vec[0] < 0 and normal_vec[1] == 0:
        return 270
    # +-/-
    elif normal_vec[1] < 0:
        return 180 + np.arctan(normal_vec[0]/normal_vec[1]) / np.pi * 180.
    # -/-
    elif normal_vec[0] < 0 and normal_vec[1] >= 0:
        return 360 + np.arctan(normal_vec[0]/normal_vec[1]) / np.pi * 180.
#    elif normal_vec[1] == 0:
#        return 90


def check_point_sets(data):
    """Checks if point sets in KmlPoints object contain at least three data points for
    fitting a plane. If not it removes the unsuitable point sets."""
    for i, ps in enumerate(data.point_sets):
        if len(ps.points) < 3:
            data.point_sets.remove(ps)
            print("Removed point set #" + str(i))


def extract_xyz(k):
    x = []
    y = []
    z = []

    for i, ps in enumerate(k.point_sets):
        for j, p in enumerate(ps.points):
            x.append(p.x)
            y.append(p.y)
            z.append(p.z)

    return np.array([x, y, z])


def points_to_gempy_interf(ks_coords, formations, series="Default series", debug=False):
    """Converts KmlPoints coordinates into GemPy interfaces dataframe.

    Args:
        ks_coords: list/array of point set arrays [[x,y,z], [x,y,z]]
        formations: list of fromation names [str, str, ...]
        series (str, optional): Set the series for the given point sets.
        debug (bool, optional): Toggles verbosity.

    Returns:
        GemPy interfaces dataframe with columns ['X', 'Y', 'Z', 'formation', 'series'].
    """
    interfaces = pn.DataFrame(columns=['X', 'Y', 'Z', 'formation', 'series'])

    for i, k in enumerate(ks_coords):
        temp = pn.DataFrame(columns=['X', 'Y', 'Z', 'formation', 'series'])
        if debug:
            print(i)

        temp["X"] = k[0]
        temp["Y"] = k[1]
        temp["Z"] = k[2]
        temp["formation"] = formations[i]
        temp["series"] = series

        interfaces = interfaces.append(temp, ignore_index=True)

    return interfaces


def dips_to_gempy_fol(dips, dip_dirs, xs, ys, zs, formation, series="Default series"):
    """

    Args:
        dips:
        dip_dirs:
        xs:
        ys:
        zs:
        formation:
        series:

    Returns:

    """
    foliations = pn.DataFrame(columns=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity', 'formation', 'series'])

    foliations["X"] = xs
    foliations["Y"] = ys
    foliations["Z"] = zs
    foliations["dip"] = dips
    foliations["azimuth"] = dip_dirs
    foliations["polarity"] = 1
    foliations["formation"] = formation
    foliations["series"] = series

    return foliations

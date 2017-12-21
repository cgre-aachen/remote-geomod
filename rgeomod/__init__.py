"""

"""

import pandas as pn
import os
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


def read_kml_files(folder_path, verbose=True):
    """Reads in all .kml files from given folder, creating a KmlPoints instance for each
    file found.
    Filename convention:
        '01_formationname_dips.kml' for dip picks.
        '04_formationname_interf.kml' for

    Args:
        folder_path (str): Relative path to the folder containing the picked points as .kml files (Google Earth).
        verbose (bool): Toggles verbosity.

    Returns:
        (list):
        (list):
        (np.array, boolean):
    """

    ks = []
    ks_names = []
    ks_bool = []

    for i, fn in enumerate(os.listdir(folder_path)):
        if ".kml" in fn:
            ks.append(kml_to_plane.KmlPoints(filename=folder_path + fn, debug=verbose))
            if verbose:
                print(fn)


            # auto check if some set contains less than 3 points and throw them out
            check_point_sets(ks[-1])
            if verbose:
                print("\n")
            # append names
            ks_names.append(fn[:-4])
            if "dips" in fn or "Dips" in fn:
                ks_bool.append(True)
            else:
                ks_bool.append(False)

    return ks, ks_names, np.array(ks_bool).astype(bool)


def fetch_z_from_dtm(ks, dtm_path):
    for k in ks:
        for ps in k.point_sets:
            try:
                ps.get_z_values_from_geotiff(dtm_path)
            except IndexError:
                print("Point outside geotiff, drop")
                k.point_sets.remove(ps)
                continue

            # convert LatLon coordinates to UTM
            ps.latlong_to_utm()
            # Fit plane to point set
            ps.plane_fit()


def calc_dips_from_points(ks, ks_bool):
    dips = []
    dip_dirs = []
    dip_xs = []
    dip_ys = []
    dip_zs = []

    for k in np.array(ks)[ks_bool]:
        for ps in k.point_sets:
            # determine dip angle from normal vector of plane
            dips.append(dip(ps.normal))
            # get dip direction from normal vector
            dip_dirs.append(dip_dir(ps.normal))
            # get centroid coordinates
            dip_xs.append(ps.ctr.x)
            dip_ys.append(ps.ctr.y)
            dip_zs.append(ps.ctr.z)

    return dips, dip_dirs, dip_xs, dip_ys, dip_zs


def convert_to_df(ks, ks_names, ks_bool):
    """

    Args:
        ks:
        ks_names:
        ks_bool:

    Returns:

    """
    # interfaces
    # ----------
    ks_coords = []
    for k in ks:
        ks_coords.append(extract_xyz(k))

    ks_coords_interf = []
    ks_names_interf = []
    for i, k in enumerate(ks_coords):
        if not ks_bool[i]:
            ks_coords_interf.append(k)
            ks_names_interf.append(ks_names[i])

    interfaces = points_to_gempy_interf(ks_coords_interf, ks_names_interf)

    # foliations
    # ----------
    dips, dip_dirs, dip_xs, dip_ys, dip_zs = calc_dips_from_points(ks, ks_bool)
    foliations = dips_to_gempy_fol(dips, dip_dirs, dip_xs, dip_ys, dip_zs, ks_names[0])

    return interfaces, foliations
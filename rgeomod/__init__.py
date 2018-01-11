"""

"""

import pandas as pn
import os
from .kml_to_plane import *
from .struct_geo import *
from tqdm import tqdm
from time import sleep
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import gdal

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


def points_to_gempy_interf(ks_coords, formations, filenames, series="Default series", debug=False):
    """Converts KmlPoints coordinates into GemPy interfaces dataframe.

    Args:
        ks_coords: list/array of point set arrays [[x,y,z], [x,y,z]]
        formations: list of fromation names [str, str, ...]
        series (str, optional): Set the series for the given point sets.
        debug (bool, optional): Toggles verbosity.

    Returns:
        GemPy interfaces dataframe with columns ['X', 'Y', 'Z', 'formation', 'series', "formation number"].
    """
    interfaces = pn.DataFrame(columns=['X', 'Y', 'Z', 'formation', 'series', "formation number"])

    for i, k in enumerate(ks_coords):
        temp = pn.DataFrame(columns=['X', 'Y', 'Z', 'formation', 'series', "formation number"])
        if debug:
            print(i)

        temp["X"] = k[0]
        temp["Y"] = k[1]
        temp["Z"] = k[2]
        temp["formation"] = formations[i]
        temp["series"] = series
        temp["formation number"] = int(filenames[i].split("_")[0])

        interfaces = interfaces.append(temp, ignore_index=True)

    return interfaces


def dips_to_gempy_fol(dips, dip_dirs, xs, ys, zs, formation, formation_number, series="Default series"):
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
    foliations["formation number"] = formation_number

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
    filenames = []

    for i, fn in enumerate(os.listdir(folder_path)):
        if ".kml" in fn:
            ks.append(kml_to_plane.KmlPoints(filename=folder_path + fn, debug=verbose))
            if verbose:
                print(fn)

            filenames.append(fn)

            # auto check if some set contains less than 3 points and throw them out
            check_point_sets(ks[-1])
            if verbose:
                print("\n")


            ks_names.append(fn.split("_")[1])  # append formation name

            if "dips" in fn or "Dips" in fn or "foliation" in fn:
                ks_bool.append(True)
            else:
                ks_bool.append(False)

    return ks, ks_names, np.array(ks_bool).astype(bool), filenames


def get_elevation_from_dtm(ks, dtm_path, verbose=True):
    for k in tqdm(ks, desc="Extracting elevation data"):
        sleep(0.3)
        for ps in k.point_sets:
            try:
                ps.get_z_values_from_geotiff(dtm_path)
            except IndexError:
                if verbose:
                    print("Point outside geotiff, drop")
                k.point_sets.remove(ps)
                continue

    if verbose:
        print("Elevation data successfully extracted from DTM.")


def fit_planes_to_points(ks, verbose=True):
    for k in tqdm(ks, desc="Fitting planes to point sets"):
        sleep(0.3)
        for ps in k.point_sets:
            # convert LatLon coordinates to UTM
            ps.latlong_to_utm()
            # Fit plane to point set
            ps.plane_fit()

    if verbose:
        print("Planes successfully fit to point sets.")


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


def convert_to_df(ks, ks_names, filenames, ks_bool):
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
    ks_filenames = []

    for i, k in enumerate(ks_coords):
        if not ks_bool[i]:
            ks_coords_interf.append(k)
            ks_names_interf.append(ks_names[i])
            ks_filenames.append(filenames[i])


    interfaces = points_to_gempy_interf(ks_coords_interf, ks_names_interf, ks_filenames)

    # foliations
    # ----------
    dips, dip_dirs, dip_xs, dip_ys, dip_zs = calc_dips_from_points(ks, ks_bool)
    foliations = dips_to_gempy_fol(dips, dip_dirs, dip_xs, dip_ys, dip_zs, ks_names[0], 1)
    # TODO: flexible formation number assignment for foliations

    return interfaces, foliations


def convert_dtm_to_gempy_grid(raster, dtm):
    # here are the raster dimensions:
    # raster.RasterXSize, raster.RasterYSize
    geoinformation = raster.GetGeoTransform()
    # get DTM corners:
    dtm_E_min = geoinformation[0]
    dtm_E_max = geoinformation[0] + geoinformation[1] * raster.RasterXSize
    dtm_N_min = geoinformation[3] + geoinformation[5] * raster.RasterYSize
    dtm_N_max = geoinformation[3]
    # dtm_E_min, dtm_E_max, dtm_N_min, dtm_N_max

    # define range for x, y - values
    X_range = np.arange(dtm_E_min, dtm_E_max, geoinformation[1])
    Y_range = np.arange(dtm_N_min, dtm_N_max, np.abs(geoinformation[5]))
    XX, YY = np.meshgrid(X_range, Y_range, indexing="ij")

    # Create list of input points for interpolation with gempy:
    return np.array(list(zip(XX.ravel(), YY.ravel(), dtm[::-1,:].T.ravel())))


def export_geotiff(path, geo_map, cmap, geotiff_filepath):
    """

    Args:
        path (str): Filepath for the exported geotiff, must end in .tif
        geo_map (np.ndarray): 2-D array containing the geological map
        cmap (matplotlib colormap): The colormap to be used for the export
        geotiff_filepath (str): Filepath of the template geotiff

    Returns:
        Saves the geological map as a geotiff to the given path.
    """

    # **********************************************************************
    geo_map_rgb = cmap(geo_map)  # r,g,b,alpha
    # **********************************************************************
    # gdal.UseExceptions()

    ds = gdal.Open(geotiff_filepath)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    [cols, rows] = arr.shape

    outFileName = path
    driver = gdal.GetDriverByName("GTiff")
    options = ['PROFILE=GeoTiff', 'PHOTOMETRIC=RGB', 'COMPRESS=JPEG']
    outdata = driver.Create(outFileName, rows, cols, 3, gdal.GDT_Byte, options=options)

    outdata.SetGeoTransform(ds.GetGeoTransform())  # sets same geotransform as input
    outdata.SetProjection(ds.GetProjection())  # sets same projection as input
    outdata.GetRasterBand(1).WriteArray(geo_map_rgb[:, ::-1, 0].T * 256)
    outdata.GetRasterBand(2).WriteArray(geo_map_rgb[:, ::-1, 1].T * 256)
    outdata.GetRasterBand(3).WriteArray(geo_map_rgb[:, ::-1, 2].T * 256)
    outdata.GetRasterBand(1).SetColorInterpretation(gdal.GCI_RedBand)
    outdata.GetRasterBand(2).SetColorInterpretation(gdal.GCI_GreenBand)
    outdata.GetRasterBand(3).SetColorInterpretation(gdal.GCI_BlueBand)
    # outdata.GetRasterBand(4).SetColorInterpretation(gdal.GCI_AlphaBand)  # alpha band

    # outdata.GetRasterBand(1).SetNoDataValue(999)##if you want these values transparent
    outdata.FlushCache()  # saves to disk
    outdata = None  # closes file (important)
    band = None
    ds = None

    print("Successfully exported geological map to "+path)


def calculate_gradient(foliations):
    """
    Calculate the gradient vector of module 1 given dip and azimuth to be able to plot the foliations

    Attributes:
        foliations: extra columns with components xyz of the unity vector.
    """

    foliations['G_x'] = np.sin(np.deg2rad(foliations["dip"].astype('float'))) * \
                             np.sin(np.deg2rad(foliations["azimuth"].astype('float'))) * \
                             foliations["polarity"].astype('float')
    foliations['G_y'] = np.sin(np.deg2rad(foliations["dip"].astype('float'))) * \
                             np.cos(np.deg2rad(foliations["azimuth"].astype('float'))) *\
                             foliations["polarity"].astype('float')
    foliations['G_z'] = np.cos(np.deg2rad(foliations["dip"].astype('float'))) *\
                             foliations["polarity"].astype('float')


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


def export_point(kmlpoints_list, template_fp, placemark_template_fp, filepath, filename):
    with open(placemark_template_fp) as f:
        placemark_template = f.readlines()

    for i, kmlp in enumerate(kmlpoints_list):
        with open(template_fp) as f:
            template = f.readlines()

        for ps in kmlp.point_sets:

            for p in ps.points:
                p.utm_to_latlong()

                # write coordinates
                placemark_template[5 - 1] = "<longitude>" + str(p.x) + "</longitude>\n"
                placemark_template[6 - 1] = "<latitude>" + str(p.y) + "</latitude>\n"
                placemark_template[7 - 1] = "<altitude>" + str(p.z) + "</altitude>\n"

                placemark_template[14 - 1] = "<coordinates>" + str(p.x) + "," + str(p.y) + "," + str(
                    p.z) + "</coordinates>\n"

                # append placemark to template
                template[-3:-3] = placemark_template

        with open(filepath+str(i)+"_"+filename+".kml", 'w') as file:
            file.write("".join(template))
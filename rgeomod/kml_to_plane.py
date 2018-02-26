"""Methods to fit plane to points picked in GoogleEarth and stored in kml files"""

try:
    from osgeo import ogr, osr
    import gdal
except ImportError:
    print("Geopgraphic libraries (osgeo, gdal) not (correctly) installed")
    print("Continuing... but some functionality may not work!")

import numpy as np


# create a class for more general calculation
class Point(object):
    def __init__(self, **kwds):
        """3-D point in space

        **Optional Keywords**:
        - 'x' = float : x or longitude
        - 'y' = float : y or latitude
        - 'z' = float : z or altitude
        - 'type' = 'utm', 'latlong', 'nongeo' # use nongeo for non-geographic projection
        - 'zone' = int: utm zone (needs to be defined for type=utm!)
        """
        # if hasattr(kwds, 'x'):
        self.x = kwds['x']
        self.y = kwds['y']
        self.type = kwds['type']
        if 'z' in kwds:
            self.z = kwds['z']
        # self.z = z
        if 'zone' in kwds:
            self.zone = kwds['zone']
        if 'type' in kwds and kwds['type'] == 'utm' and not 'zone' in kwds:
            raise AttributeError("Please provide utm zone")
        # for the case that self.type == 'latlong': determine UTM zone:
        if self.type == 'latlong':
            self.zone = int(np.floor(np.mod((self.x + 180) / 6, 60)) + 1)

    def __repr__(self):
        if hasattr(self, 'z'):
            return "p(%f, %f, %f) in %s" % (self.x, self.y, self.z, self.type)
        else:
            return "p(%f, %f) in %s" % (self.x, self.y, self.type)

    def latlong_to_utm(self):
        """Convert point from lat long to utm for given zone"""
        wgs = osr.SpatialReference()
        wgs.ImportFromEPSG(4326)
        if self.zone == 40:
            utm = osr.SpatialReference()
            utm.ImportFromEPSG(32640)
        else:
            raise AttributeError("Sorry, zone %d not yet implemented\
             (to fix: check EPSG code on http://spatialreference.org/ref/epsg/ and include in code!)" % self.zone)
        ct = osr.CoordinateTransformation(wgs, utm)
        self.x, self.y = ct.TransformPoint(self.x, self.y)[:2]
        self.type = 'utm'

    def utm_to_latlong(self):
        """Convert point from utm to lat long for given zone"""
        wgs = osr.SpatialReference()
        wgs.ImportFromEPSG(4326)
        if self.zone == 40:
            utm = osr.SpatialReference()
            utm.ImportFromEPSG(32640)
        else:
            raise AttributeError("Sorry, zone %d not yet implemented (check EPSG code and include in code!)" % self.zone)
        ct = osr.CoordinateTransformation(utm, wgs)
        self.x, self.y = ct.TransformPoint(self.x, self.y)[:2]
        self.type = 'latlong'

    def get_z_value_from_geotiff(self, filename):
        """Load GeoTiff and extract pixel value at position of point"""
        if self.type == 'utm':
            self.utm_to_latlong(self.zone)

        l = looker(filename)
        self.z = l.lookup(self.x, self.y)


class PointSet(object):

    def __init__(self, **kwds):
        """Point set as a collection of points (picks on one line)

        **Optional keywords**:
            - *type* = 'utm', 'latlong': coordinate type (default: 'latlong')
        """
        self.points = []
        self.type = kwds.get('type', 'latlong')

    def __repr__(self):
        """Print out information about point set"""
        str = "Point set with %d points" % len(self.points)
        str += "; " + self.type
        if hasattr(self, 'ctr'):
            str+="; Centroid: at (%.2f, %.2f, %.2f)" % (self.ctr.x, self.ctr.y, self.ctr.z)

        if hasattr(self, 'dip'):
            str+="; Orientation: (%03d/%02d)" % (self.dip_direction, self.dip)
        return str

    def add_point(self, point):
        self.points.append(point)

    def latlong_to_utm(self):
        """Convert all points from lat long to utm"""
        if self.type == 'latlong': # else not required...
            for point in self.points:
                point.latlong_to_utm()
            self.type = 'utm'
        # convert plane centroid, if already calculated:
        if hasattr(self, 'ctr'):
            self.ctr.latlong_to_utm()

    def utm_to_latlong(self):
        """Convert all points from utm to lat long"""
        if self.type == 'utm': # else not required...
            for point in self.points:
                point.utm_to_latlong()
            self.type = 'latlong'
        # convert plane centroid, if already calculated:
        if hasattr(self, 'ctr'):
            self.ctr.utm_to_latlong()

    def get_z_values_from_geotiff(self, filename):
        """Open GeoTiff file and get z-value for all points in set"""

        # check if points in latlong, else: convert
        if self.type == 'utm':
            self.utm_to_latlong(self.zone)

        # initialise lookup for entire point set
        l = looker(filename)

        for point in self.points:
            point.z = l.lookup(point.x, point.y)

    def plane_fit(self):
        """Fit plane to points in PointSet

        Fit an d-dimensional plane to the points in a point set
        Return a point, p, on the plane (the point-cloud centroid),
        and the normal, n.

        adjusted from: http://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
        """
        import numpy as np

        if self.type == 'latlong':
            self.latlong_to_utm()

        points = np.empty((3, len(self.points)))
        for i, point in enumerate(self.points):
            points[0, i] = point.x
            points[1, i] = point.y
            points[2, i] = point.z

        from numpy.linalg import svd
        points = np.reshape(points, (np.shape(points)[0], -1))  # Collapse trialing dimensions
        assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1],
                                                                                                       points.shape[0])
        ctr = points.mean(axis=1)
        x = points - ctr[:, np.newaxis]
        M = np.dot(x, x.T)  # Could also use np.cov(x) here.

        self.ctr = Point(x=ctr[0], y=ctr[1], z=ctr[2], type='utm', zone=self.points[0].zone)
        self.normal = svd(M)[0][:, -1]
        # return ctr, svd(M)[0][:, -1]
        if self.normal[2] < 0:
            self.normal = - self.normal

    def get_orientation(self):
        """Get orientation (dip_direction, dip) for points in all point set"""
        if "normal" not in dir(self):
            self.plane_fit()

        # calculate dip
        self.dip = np.arccos(self.normal[2]) / np.pi * 180.

        # calculate dip direction
        # +/+
        if self.normal[0] >= 0 and self.normal[1] > 0:
            self.dip_direction = np.arctan(self.normal[0]/self.normal[1]) / np.pi * 180.
        # border cases where arctan not defined:
        elif self.normal[0] > 0 and self.normal[1] == 0:
            self.dip_direction = 90
        elif self.normal[0] < 0 and self.normal[1] == 0:
            self.dip_direction = 270
        # +-/-
        elif self.normal[1] < 0:
            self.dip_direction = 180 + np.arctan(self.normal[0]/self.normal[1]) / np.pi * 180.
        # -/-
        elif self.normal[0] < 0 and self.normal[1] >= 0:
            self.dip_direction = 360 + np.arctan(self.normal[0]/self.normal[1]) / np.pi * 180.
#    elif normal_vec[1] == 0:
#        return 90

    def stereonet(self):
        """Create stereonet plot of plane pole and half circle for this point set"""
        import matplotlib.pyplot as plt
        import mplstereonet

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='stereonet')

        # ax.plane(dip_dirs, dips, 'g-', linewidth=0.5)
        ax.pole(self.dip_direction - 90, self.dip, 'gs', markersize=4)
        ax.plane(self.dip_direction - 90, self.dip, 'g', markersize=4)

        # ax.rake(strike, dip, -25)
        ax.grid()

    def minmax(self):
        """Get minimum and maximum values of points in point set (e.g. to determine surrounding box)"""
        point_array = np.empty((len(self.points), 2))

        for i, p in enumerate(self.points):
            point_array[i] = (p.x, p.y)

        self.min = np.min(point_array, axis=0)
        self.max = np.max(point_array, axis=0)


#
#   NOTE: code for plane is somewhat deprecated - plane fit is based on z=f(x,y) assumption
#
class Plane(object):

    def __init__(self):
        """Plane as result of point fit"""
        self.points = []

    def add_point(self, Point):
        """Add point to plane"""
        self.points.append(Point)

    def update_coefficients(self):
        """Determine plane coefficients"""
        # assemble matrices:
        self.B = np.array([[1, p.x, p.y] for p in self.points])
        self.z = np.array([p.z for p in self.points])
        self.a = np.dot(np.linalg.pinv(self.B), self.z)

    def update_orientation(self):
        """Update geological orientation values"""
        self.update_coefficients()
        # calculate dip and azimuth:
        maxslope = np.sqrt(self.a[1] ** 2 + self.a[2] ** 2)
        self.dip = np.arctan(maxslope) / np.pi * 180.
        self.dip_direction = np.abs(np.arctan2(-self.a[1], -self.a[2])) / np.pi * 180.

    def get_orientation(self):
        """Get geological orientation values of plane"""
        self.update_orientation()
        return self.dip, self.dip_direction

    def get_misfit(self):
        """Get sum of squared misfit"""
        self.misfit = np.sum((np.dot(self.B, self.a) - self.z) ** 2)
        return self.misfit

    def __repr__(self):
        """Print information about plane"""
        str = "Plane fit through %d points\n" % len(self.points)
        if hasattr(self, 'dip'):
            str += "Orientation: (%06.2f, %05.2f)\n" % (self.dip_direction, self.dip)
        if hasattr(self, 'misfit') and len(self.points) > 3:
            # show misfit only if more than three points
            str += ("Misfit (L2): %.5f\n" % self.misfit)
        return str



def planeFit(PointSet):
    """
    p, n = planeFit(PointSet)

    Fit an d-dimensional plane to the points in a point set
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.

    adjusted from: http://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    """
    import numpy as np

    points = np.empty((3, len(PointSet.points)))
    for i, point in enumerate(PointSet.points):
        points[0, i] = point.x
        points[1, i] = point.y
        points[2, i] = point.z

    import numpy as np
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1))  # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1],
                                                                                                   points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    M = np.dot(x, x.T)  # Could also use np.cov(x) here.
    return ctr, svd(M)[0][:, -1]


class KmlPoints(object):

    def __init__(self, **kwds):
        """Get point sets from KML file

        **Optional keywords**:

        - *filename* = string: filename of kml file
        - *debug* = bool: provide debug output (Default: false)
        - *auto_remove* = bool: automatically remove unsuitable points (e.g. outside Geotiffs)
            and point sets (e.g. too few points, too close on a line)
        - 'type' = 'utm', 'latlong' : coordinate system of points (default: latlong)
        """
        self.debug = kwds.get("debug", False)
        self.auto_remove = kwds.get("auto_remove", True)
        self.type = kwds.get("type", 'latlong')
        self.geotiffs = []
        self.points = []
        self.point_sets = []
        # if kwds.has_key('filename'):
        if 'filename' in kwds:
            if self.debug:
                print("read kml")
            self.read_kml(kwds['filename'])

    def read_kml(self, filename):
        """Read kml file and extract points"""

        ds = ogr.Open(filename)
        point_sets = []

        for lyr in ds:
            for j, feat in enumerate(lyr):
                geom = feat.GetGeometryRef()
                ps = PointSet()
                if geom != None:
                    for i in range(0, geom.GetPointCount()):
                        # print (geom.GetPoint(i))
                        point = Point(x = geom.GetPoint(i)[0],
                                      y = geom.GetPoint(i)[1],
                                      type = 'latlong')
                        ps.add_point(point)
                        # points.append([geom.GetPoint(i)[0], geom.GetPoint(i)[1], j])

                self.point_sets.append(ps)

        if self.debug:
            print("%d point sets added" % len(self.point_sets))

    def test_point_sets(self):
        """Test if point sets contain at least three points; if not: remove"""
        # test if all point sets have at least three points:
        for ps in self.point_sets:
            if len(ps.points) < 3:
                self.point_sets.remove(ps)
                if self.debug:
                    print("Removed point set")

        if self.debug:
            print("%d point sets remaining" % len(self.point_sets))

    def determine_z_values(self):
        """Determine z values for all points in point sets

        Approach: test all geotiffs in given order, stored in self.geotiffs list
        """
        if len(self.geotiffs) == 0:
            raise AttributeError("Please define geotiffs first (self.add_geotiff())")

        # check that coordinates are in latlong, if not: convert
        if self.type == 'utm':
            self.utm_to_latlong()

        for ps in self.point_sets:
            fail = True
            for geotiff in self.geotiffs:
                try:
                    ps.get_z_values_from_geotiff(geotiff)
                except IndexError:
                    continue
                fail = False

            # if point can not be detected: remove (default) or raise error
            # if self.auto_remove = False

            if fail:
                if self.auto_remove:
                    if self.debug:
                        print("Point outside geotiff, drop")
                    self.point_sets.remove(ps)
                else:
                    raise IndexError("Point outside of defined geotiffs!\nPlease define\
                                     suitable geotiff or remove point (set self.auto_remove = True)")

    def fit_plane_to_all_sets(self):
        """Fit plane to all point sets

        Results are stored in point set object (self.ctr, self.normal)
        """
        if self.type == 'latlong':
            self.latlong_to_utm()

        for ps in self.point_sets:
            ps.plane_fit()
            ps.get_orientation()

    def stereonet(self):
        """Create stereonet plot of all plane pole and half circle for all planes"""
        import matplotlib.pyplot as plt
        import mplstereonet

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='stereonet')

        # ax.plane(dip_dirs, dips, 'g-', linewidth=0.5)
        for ps in self.point_sets:
            ax.pole(ps.dip_direction - 90, ps.dip, 'gs', markersize=4)
            ax.plane(ps.dip_direction - 90, ps.dip, 'g', markersize=4)

        # ax.rake(strike, dip, -25)
        ax.grid()

    def export_for_geomodeller(self, filename='points.csv', formname="formation", data_type='ori', **kwds):
        """Export points to geomodeller

        **Optional Keywords**:

            - *filename* = string: file name of output csv file (default: points.csv)
            - *formname* = string: name of geological formation
            - *data_type* = 'planes' / 'ori' : type of picks (plane tops, orientation measurements); default: 'ori'
        """
        # initialise file
        filename = filename
        formname = formname
        data_type = data_type

        f = open(filename, 'w')

        if data_type == 'ori':
            f.write("X,Y,Z,azimuth,dip,polarity,formation\n")

            for ps in self.point_sets:
                f.write("%.1f,%.1f,%.1f,%.2f,%.2f,1,%s\n" % \
                        (ps.ctr.x, ps.ctr.y, ps.ctr.z, ps.dip_direction, ps.dip, formname))
            f.close()

        elif data_type == 'planes':
            f.write("X,Y,Z,formation\n")
            for ps in self.point_sets:
                for p in ps.points:
                    f.write("%.1f,%.1f,%.1f,%s\n" % \
                        (p.x, p.y, p.z, formname))
            f.close()
        # else:
         #    raise AttributeError("data type not known")


    def add_geotiff(self, geotiff):
        """Add geotiff to list of geotiffs (self.geotiffs)

        **Arguments**:
            - *geotiff* = filename : filename (with complete path) to geotiff
        """
        if self.debug:
            print("Note: for efficiency reasons, add the most important geotiff first!")
        self.geotiffs.append(geotiff)

    def export_ori_for_centroids(self):
        """Export orientation information for centroids of point set

        Format in simple csv file, usable with Geomodeller
        """
        pass

    def export_ori_all_points(self):
        """Exort orientation information for all points in point set

        Format in simple csv file, usable with Geomodeller
        """
        pass

    def latlong_to_utm(self):
        """Convert all points from lat long to utm"""
        if self.type == 'latlong':  # else not required...
            if self.debug:
                print("Convert Lat/Long to UTM")
            for ps in self.point_sets:
                ps.latlong_to_utm()
            self.type = 'utm'

    def utm_to_latlong(self):
        """Convert all points from utm to lat long"""
        if self.type == 'utm':  # else not required...
            if self.debug:
                print("Convert UTM to Lat/Long")
            for ps in self.point_sets:
                ps.utm_to_latlong()
            self.type = 'latlong'


class looker(object):
    """let you look up pixel value

    Credits to entry on stackoverflow:
    http://stackoverflow.com/questions/13439357/extract-point-from-raster-in-gdal
    """

    def __init__(self, tifname='test.tif'):
        """Give name of tif file (or other raster data?)"""

        # open the raster and its spatial reference
        self.ds = gdal.Open(tifname)
        srRaster = osr.SpatialReference(self.ds.GetProjection())

        # get the WGS84 spatial reference
        srPoint = osr.SpatialReference()
        srPoint.ImportFromEPSG(4326) # WGS84

        # coordinate transformation
        self.ct = osr.CoordinateTransformation(srPoint, srRaster)

        # geotranformation and its inverse
        gt = self.ds.GetGeoTransform()
        dev = (gt[1]*gt[5] - gt[2]*gt[4])
        gtinv = ( gt[0] , gt[5]/dev, -gt[2]/dev,
                gt[3], -gt[4]/dev, gt[1]/dev)
        self.gt = gt
        self.gtinv = gtinv

        # band as array
        b = self.ds.GetRasterBand(1)
        self.arr = b.ReadAsArray()

    def lookup(self, lon, lat):
        """look up value at lon, lat"""

        # get coordinate of the raster
        xgeo,ygeo,zgeo = self.ct.TransformPoint(lon, lat, 0)

        # convert it to pixel/line on band
        u = xgeo - self.gtinv[0]
        v = ygeo - self.gtinv[3]
        # FIXME this int() is probably bad idea, there should be
        # half cell size thing needed
        xpix = int(self.gtinv[1] * u + self.gtinv[2] * v)
        ylin = int(self.gtinv[4] * u + self.gtinv[5] * v)

        # look the value up
        return self.arr[ylin, xpix]


# template from kml geosymbol generator: http://csmres.jmu.edu/Geollab/Whitmeyer/web/visuals/GoogleEarth/tools/SD.html

# Set variables:
# %dip%
# %strike%
# %lat%
# %long%
# %name%
#


kml_template = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
	<name>Orientation symbols</name>
	<Style id="sn_no_icon">
		<IconStyle><Icon></Icon></IconStyle>
		<LabelStyle>
			<scale>1.0</scale>
		</LabelStyle>
	</Style>
	<Style id="sn_shaded_dot">
		<IconStyle>
			<color>03ffffff</color>
			<scale>1.2</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>
			</Icon>
		</IconStyle>
		<LabelStyle>
			<color>00ffffff</color>
		</LabelStyle>
		<BalloonStyle>
			<text>$[description]</text>
		</BalloonStyle>
	</Style>
	<Placemark>
		<name>Symbol1</name>
		<Model id="SDmodel">
			<altitudeMode>relativeToGround</altitudeMode>
			<Location>
				<longitude>%long%</longitude>
				<latitude>%lat%</latitude>
				<altitude>20</altitude>
			</Location>
			<Orientation>
				<heading>%strike%</heading>
				<tilt>0</tilt>
				<roll>-%dip%</roll>
			</Orientation>
			<Scale>
				<x>40</x>
				<y>50</y>
				<z>50</z>
			</Scale>
			<Link>
				<href>http://csmres.jmu.edu/Geollab/Whitmeyer/web/visuals/GoogleEarth/tools/SDblack.dae</href>
			</Link>
		</Model>
	</Placemark>
	<Placemark>
		<name>test1</name>
		<description>Unit: test
Strike: %strike%
Dip: %dip$

Notes:
</description>
		<styleUrl>#sn_shaded_dot</styleUrl>
		<Point>
			<coordinates>%long%,%lat%,0</coordinates>
		</Point>
	</Placemark>
	<Placemark>
		<name>%dip$</name>
		<styleUrl>#sn_no_icon</styleUrl>
		<Point>
			<coordinates>%long%,%lat%,0</coordinates>
		</Point>
	</Placemark>
</Document>
</kml>"""
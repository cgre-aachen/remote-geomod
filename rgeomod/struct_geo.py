"""Functions and classes for basic structural geology and mapping applications"""

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def calculate_dip_direction_dip(p_0, p_1, p_2):
    """Calculate dip direction and dip of a plane

    **Arguments**:

    p_0 = (x_0, y_0, z_0) : coordinates of first point
    p_1 = (x_1, y_1, z_1) : coordinates of second point
    p_2 = (x_2, y_2, z_2) : coordinates of third point

     **Returns**:

    dip_direction, dip : tuple of two numbers
    """
    # YOUR CODE HERE
    # create an array containing all points and extract x,y and z
    points = np.column_stack((np.array(p_0), np.array(p_1), np.array(p_2)))

    B = np.ones((3, 3))
    B[:, 1] = points[0, :]
    B[:, 2] = points[1, :]

    z = points[2, :]

    # solving the linear interpolation for a
    (a_0, a_1, a_2) = np.linalg.solve(B, z)

    maxslope = np.sqrt(a_1 ** 2 + a_2 ** 2)
    dip = np.arctan(maxslope) / np.pi * 180.
    dip_direction = np.abs(np.arctan2(-a_1, -a_2)) / np.pi * 180.

    return (dip, dip_direction)


def plot_three_point(p_0, p_1, p_2):
    """Interpolate plane between three points and generate plot"""

    import matplotlib.patches as mpatches

    # extracting x,y,z values again
    points = np.column_stack((np.array(p_0), np.array(p_1), np.array(p_2)))
    x_vals = points[0, :]
    y_vals = points[1, :]
    z_vals = points[2, :]

    # repeating the linear interpolation
    B = np.ones((3, 3))
    B[:, 1] = x_vals
    B[:, 2] = y_vals

    a = np.linalg.solve(B, z_vals)

    # creating a meshgrid of given x,y values and calculating Z_grid
    x_grid = np.arange(-2, 2.5, 0.5)
    y_grid = np.arange(-2, 2.5, 0.5)

    (X_grid, Y_grid) = np.meshgrid(x_grid, y_grid)

    Z_grid = a[0] + X_grid * a[1] + Y_grid * a[2]

    # obtaining the dip and dip direction for annotation in the plot
    (dip, dip_dir) = calculate_dip_direction_dip(p_0, p_1, p_2)

    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    plot_points, = ax.plot(x_vals, y_vals, z_vals, 'o', color='#19647E', label='Original points')
    ax.plot_surface(X_grid, Y_grid, Z_grid, color='#119DA4', shade=False)
    # labeling of the surface did not work, creating a patch for the legend (proxy artist)
    surface_patch = mpatches.Patch(color='#119DA4', label='Interpolated plane')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Determination of a plane with three points")
    ax.view_init(20, 100)
    ax.annotate(xy=(0, 0),
                xytext=(0.01, -0.01),
                bbox=dict(boxstyle='round,pad=0.4', fc='#FFC857'),
                s='Dip direction / dip: ( %03d | %d )' % (round(dip_dir), round(dip)),
                color='#000000',
                xycoords='data'
                )
    ax.legend(handles=[plot_points, surface_patch])


def plane_fit(points):
    """Fit plane to points in PointSet

    Fit an d-dimensional plane to the points in a point set
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.

    adjusted from: http://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    """
    import numpy as np

    # points = np.empty((3, len(points)))

    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1))  # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1],
                                                                                                   points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    M = np.dot(x, x.T)  # Could also use np.cov(x) here.

    # self.ctr = Point(x=ctr[0], y=ctr[1], z=ctr[2], type='utm', zone=self.points[0].zone)
    normal = svd(M)[0][:, -1]
    # return ctr, svd(M)[0][:, -1]
    if normal[2] < 0:
        normal = - normal

    a_1, a_2, a_3, = normal

    maxslope = np.sqrt(a_1 ** 2 + a_2 ** 2)
    dip = np.arctan(maxslope) / np.pi * 180.
    dip_direction = np.abs(np.arctan2(-a_1, -a_2)) / np.pi * 180.

    return (dip, dip_direction)

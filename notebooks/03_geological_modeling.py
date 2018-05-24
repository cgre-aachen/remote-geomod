
# coding: utf-8

# # 3 - 3D Modeling with GemPy

# In[1]:

from matplotlib import use
use("Agg")

import sys
import numpy as np
# These two lines are necessary only if gempy is not installed
sys.path.append("../../gempy/")
sys.path.append("../gempy/")

# Importing gempy
import gempy as gp
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


geo_data=gp.create_data(extent=[612000, 622000, 2472000, 2480000, -1000, 1000], 
                        resolution=[50, 50, 50], 
                        path_f = "../data/gempy_foliations.csv",
                        path_i = "../data/gempy_interfaces.csv")


# In[3]:


gp.set_series(geo_data, {"Default series": np.unique(geo_data.interfaces["formation"].values)},
             order_formations = np.unique(geo_data.interfaces["formation"].values))

gp.set_order_formations(geo_data, np.unique(geo_data.interfaces["formation"].values))


# ## 3.2 - Data visualization

# In[4]:


gp.plot_data(geo_data, direction="z")


# In[5]:


gp.plot_data(geo_data, direction="x")


# ## 3.3 - Computing the 3D Model

# In[ ]:


interp_data = gp.InterpolatorInput(geo_data, dtype="float32")


# In[ ]:


lith_block, fault_block = gp.compute_model(interp_data)


# ## 3.4 - Model visualization
# 
# ### 3.4.1 - 2D Sections

# In[ ]:


gp.plot_section(geo_data, lith_block[0], 25, direction='y')


# In[ ]:


gp.plot_section(geo_data, lith_block[0], 25, direction='x')


# ### 3.4.2 - Pseudo-3D surfaces

# In[ ]:


v_l, s_l = gp.get_surfaces(interp_data, potential_lith=lith_block[1], step_size=2)


# In[ ]:


fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111, projection='3d')

cs = ["blue", "red", "orange"]



for i in range(3):
    surf = ax.plot_trisurf(v_l[i][:,0], v_l[i][:,1], v_l[i][:,2], 
                           color=cs[i], linewidth=0, alpha=0.75, shade=False)


# # Exporting a geological map

# In[ ]:


import gdal
geotiff_filepath = "../data/dome_sub_sub_utm.tif"
raster = gdal.Open(geotiff_filepath)
dtm = raster.ReadAsArray()


# In[ ]:


import matplotlib.pyplot as plt
plt.imshow(dtm)
plt.show()


# In[ ]:


# here are the raster dimensions:
raster.RasterXSize, raster.RasterYSize
geoinformation = raster.GetGeoTransform()


# get DTM corners:
# 

# In[ ]:


dtm_E_min = geoinformation[0]
dtm_E_max = geoinformation[0] + geoinformation[1] * raster.RasterXSize
dtm_N_min = geoinformation[3] + geoinformation[5] * raster.RasterYSize
dtm_N_max = geoinformation[3]
dtm_E_min, dtm_E_max, dtm_N_min, dtm_N_max


# In[ ]:


# define range for x, y - values
X_range = np.arange(dtm_E_min, dtm_E_max, geoinformation[1])
Y_range = np.arange(dtm_N_min, dtm_N_max, np.abs(geoinformation[5]))
XX, YY = np.meshgrid(X_range, Y_range)


# Create list of input points for interpolation with gempy:

# In[ ]:


points = np.array(list(zip(XX.ravel(), YY.ravel(), dtm.ravel())))


# In[ ]:


points


# Build basic gempy using _data_ extent (_resolution_ does not matter, as far as I can tell):

# In[ ]:


geo_data=gp.create_data(extent=[612000, 622000, 2472000, 2480000, -1000, 1000],  
                        resolution=[1, 1, 1], 
                        path_f = "../data/gempy_foliations.csv",
                        path_i = "../data/gempy_interfaces.csv")

gp.set_series(geo_data, {"Default series": np.unique(geo_data.interfaces["formation"].values)},
             order_formations = np.unique(geo_data.interfaces["formation"].values))

gp.set_order_formations(geo_data, np.unique(geo_data.interfaces["formation"].values))



# Now here the "trick": replace grid points with DTM grid points:

# In[ ]:


geo_data.grid.grid = points


# Perform the "usual" interpolation step:

# In[ ]:


interp_data = gp.InterpolatorInput(geo_data, dtype="float32")
lith_block, fault_block = gp.compute_model(interp_data)


# And here: **the geological map**:

# In[ ]:


geo_map = lith_block[0].copy().reshape((271,339))
# adjust scale
geo_map += 1
# adjust misfits:
geo_map[np.where(geo_map==5)] = 1
geo_map[np.where(geo_map==6)] = 2
geo_map[np.where(geo_map==7)] = 3
geo_map[np.where(geo_map==8)] = 4
geo_map[np.where(geo_map==1)] = 5
geo_map -= 1
# change to int for later use:
geo_map = geo_map.astype('int16')
# adjust orientation:
geo_map = geo_map[::-1,:]
plt.imshow(geo_map, cmap='viridis')
plt.colorbar()


# In[ ]:


import matplotlib.pyplot as plt
plt.imshow(dtm)
plt.show()


# In[ ]:


np.min(dtm), np.max(dtm)


# In[ ]:


file = geotiff_filepath
outFileName = "/data/geomap.tif"
# transform data
# geo_map = geo_map.astype('int16')
ds = gdal.Open(file)
band = ds.GetRasterBand(1)
arr = band.ReadAsArray()
[cols, rows] = arr.shape
# arr_min = arr.min()
# arr_max = arr.max()
# arr_mean = int(arr.mean())
# arr_out = np.where((arr < arr_mean), 10000, arr)
driver = gdal.GetDriverByName("GTiff")
# options = ['PHOTOMETRIC=RGB', 'PROFILE=GeoTIFF']
options = ['PROFILE=GeoTiff', 'PHOTOMETRIC=RGB', 'COMPRESS=JPEG']
# outdata = driver.Create(outFileName, rows, cols, 3, gdal.GDT_UInt16, options=options)
outdata = driver.Create(outFileName, rows, cols, 3, gdal.GDT_Byte, options=options)


outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
outdata.SetProjection(ds.GetProjection())##sets same projection as input
outdata.GetRasterBand(1).WriteArray(geo_map*64)
outdata.GetRasterBand(2).WriteArray(geo_map*64)
outdata.GetRasterBand(3).WriteArray(geo_map*64)
# outdata.GetRasterBand(4).WriteArray(np.ones_like(geo_map)*100)
outdata.GetRasterBand(1).SetColorInterpretation(gdal.GCI_RedBand)
outdata.GetRasterBand(2).SetColorInterpretation(gdal.GCI_GreenBand)
outdata.GetRasterBand(3).SetColorInterpretation(gdal.GCI_BlueBand)
# outdata.GetRasterBand(4).SetColorInterpretation(gdal.GCI_AlphaBand)


# outdata.GetRasterBand(1).SetNoDataValue(999)##if you want these values transparent
outdata.FlushCache() ##saves to disk!!
outdata = None
band=None
ds=None


# In[ ]:


plt.imshow(geo_map)
geo_map[np.where(geo_map == 0)] = 10
plt.colorbar()
plt.show()


# In[ ]:


arr.dtype


# In[ ]:


geo_map.dtype


# In[ ]:


ds = gdal.Open(file)


# In[ ]:


ds.GetRasterBand


# In[ ]:


np.min(geo_map), np.max(geo_map*64)


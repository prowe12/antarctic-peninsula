"""


Installing IRIS 
$ conda install -c conda-forge iris 
"""

import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

main_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/"
data_dir = main_dir + "model_performance_KGI_2017_01/models/MetUM/"
filename = "1p5km/20170131T0000Z_Escudero_1p5km_m01s16i004_air_temperature.nc"
# filename = '12km/20170131T0000Z_Escudero_12km_m01s16i004_air_temperature.nc'


# Location of Escudero station
lat = -62 - 12 / 60 - 5 / 60 / 60  # 62o 12' 05'' S
lon = -58 - 57 / 60 - 44 / 60 / 60  # 58o 57' 44'' W


# read all the variables from the 12km air temperature (on model levels) file
cubes = iris.load(data_dir + filename)

# show what was read - the netCDF files for variables on model levels include
# "surface_altitude" so that the altitude (which varies by level and grid point
# as the model levels are hybrid height coordinates) can be derived
print(cubes)

# extract the air_temperature (which could also be extracted by index, using ta = cubes[0])
ta = cubes.extract("air_temperature")[0]

# show the shape of the cube
print(ta.shape)

# which is also the shape of the cube's data (which is read lazily - so this line
# causes the data to be read)
print(ta.data.shape)

# print the cube to see its shape, coordinates and attributes
print(ta)

# print the time and model_level_number coordinates
print(ta.coord("time"))
print(ta.coord("model_level_number"))

# or access time and model_level_number coordinates using their CF coordinate types
# (see http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#coordinate-types)
#
# note the use of dim_coords=True so that only dimension coordinates are
# matched as the cube also has auxiliary and derived dimensions spanning the
# T and Z axes (shown when the cube is printed)
print(ta.coord(axis="T", dim_coords=True))
print(ta.coord(axis="Z", dim_coords=True))

# print the dimensions spanned by the altitude auxiliary coordinate
altitude = ta.coord("altitude")
for i, coord_index in enumerate(ta.coord_dims(altitude)):
    dim_coord = ta.dim_coords[coord_index]
    print(
        "Dimension {} of coordinate {} is {} (size {})".format(
            i, altitude.name(), dim_coord.name(), len(dim_coord.points)
        )
    )

# show the range of altitudes spanned by data on each model level
model_levels = ta.coord("model_level_number")
for i, level in enumerate(model_levels.points):
    level_altitudes = altitude.points[i, ...]
    print(
        "Model level {} (index {}): {} to {} {}".format(
            level,
            i,
            level_altitudes.min(),
            level_altitudes.max(),
            altitude.units,
        )
    )

# extract the temperature on the lowest level above 5000 m ASL
# (model level 38 - see above) for 2017-01-31 12Z
#
# note that a time constraint of time=iris.time.PartialDateTime(hour=12) would have
# also worked as we only have hourly data from 03Z to 00Z the following day
#
# note that the cube can be extracted by indexing, as ta_to_plot = ta[9, 37]
# but this involves knowing the order of the cube's dimensions and the relevant indices
ta_to_plot = ta.extract(
    iris.Constraint(
        model_level_number=38, time=iris.time.PartialDateTime(2017, 1, 31, 12)
    )
)

# show what we will plot
print(ta_to_plot)

# plot the extracted data on a polar stereographic grid with the domain's central
# longitude shown vertically
#
# get the coordinate reference system (CRS) for the data, and set up the CRS
# for the plot
data_cs = ta_to_plot.coord_system()
plot_cs = iris.coord_systems.Stereographic(
    central_lat=-90.0,
    central_lon=data_cs.grid_north_pole_longitude,
    true_scale_lat=-71.0,
    ellipsoid=data_cs.ellipsoid,
)

fig = plt.figure(figsize=[6, 7])

# Note: specifying a projection that is a cartopy CRS means that the Axes
# object created is actually a cartopy.mpl.geoaxes.GeoAxes object that knows
# the target projection
ax = plt.axes([0.05, 0.0, 0.9, 0.9], projection=plot_cs.as_cartopy_crs())

# iris.quickplot.contourf() is a convenience wrapper of iris.plot.contourf()
# that also gives the plot a title and a colorbar. iris.plot.contourf()
# calls matplotlib.pyplot.contourf() with a transform object that is a
# cartopy version of the cube's coordinate system, so it knows both the source
# and target coordinate systems and can project the data appropriately
#
# note that iris.plot/quickplot.contourf() take optional arguments to define
# contour levels, color map and so on (passed through to
# matplotlib.pyplot.contourf())
qplt.contourf(ta_to_plot)

# add detailed coastlines and grid lines to the plot, and show it
ax.coastlines(resolution="10m")
ax.gridlines()
plt.show()


lon = 180.0
dlon = 0.0135 - 0.0001
lat = 0.0
dlat = 0.0135 - 0.0001
this_time = iris.time.PartialDateTime(2017, 1, 31, 12)

lons = lambda cell: lon - dlon < cell < lon + dlon
lats = lambda cell: lat - dlat < cell < lat + dlat
times = iris.time.PartialDateTime(2017, 1, 31, 12)


# Get the temperature at a given lat/lon/time
time_cons = iris.Constraint(time=this_time)
lon_cons = iris.Constraint(grid_longitude=lons)
lat_cons = iris.Constraint(grid_latitude=lats)
theta = ta.extract(time_cons & lon_cons & lat_cons)


# Plot these profiles on the same set of axes. In each case we call plot
# with two arguments, the cube followed by the depth coordinate. Putting
# them in this order places the depth coordinate on the y-axis.
# The first plot is in the default axes. We'll use the same color for the
# curve and its axes/tick labels.
temperature_color = "blue"
plt.figure(figsize=(5, 6))
ax1 = plt.gca()
iplt.plot(theta, theta.coord("model_level_number"))
ax1.set_xlabel("Temperature (K)", color=temperature_color)
ax1.set_ylabel("Model level")
for ticklabel in ax1.get_xticklabels():
    ticklabel.set_color(temperature_color)


import cartopy

cartopy.crs.CRS.transform_point(lat, lon)

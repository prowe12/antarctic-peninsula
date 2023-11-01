import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

main_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/"
data_dir = main_dir + "model_performance_KGI_2017_01/models/MetUM/"
filename = "1p5km/20170131T0000Z_Escudero_1p5km_m01s16i004_air_temperature.nc"

# read air temperature from the 1.5 km run
ta = iris.load_cube(data_dir + filename, "air_temperature")

# extract data for 18Z on 31-Jan-2017
ta = ta.extract(
    iris.Constraint(time=iris.time.PartialDateTime(2017, 1, 31, 18))
)

# point of interest: 62.2S, 58.97W
lat = -62.2
lon = -58.97

# get grid coordinates corresponding to point of interest
data_cs = ta.coord_system()
latlon_cs = iris.coord_systems.GeogCS(
    semi_major_axis=data_cs.ellipsoid.semi_major_axis
)
rlonlat = data_cs.as_cartopy_crs().transform_point(
    lon, lat, latlon_cs.as_cartopy_crs()
)
rlon = rlonlat[0]
rlat = rlonlat[1]

# extract data at the point of interest
ta_at_point = ta.interpolate(
    [("grid_longitude", rlon), ("grid_latitude", rlat)],
    iris.analysis.Linear(extrapolation_mode="mask"),
)

# plot the temperature profile
qplt.plot(ta_at_point, ta_at_point.coord("altitude"))

# update the title to show point and time
title = plt.gca().title
tc = ta_at_point.coord("time")
title.set_text(
    "{} at {}\N{DEGREE SIGN}{}, {}\N{DEGREE SIGN}{} at {}".format(
        title.get_text(),
        abs(lat),
        "N" if lat >= 0 else "S",
        abs(lon),
        "E" if lon >= 0 else "W",
        tc.units.num2date(tc.points[0]).strftime(),
    )
)

# Plot temperature from radiosonde

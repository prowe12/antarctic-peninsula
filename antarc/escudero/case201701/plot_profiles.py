import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt

# My modules
from antarc.era5_read_profiles import read_era5_profiles_interp

# Directories
main_dir = "/Users/prowe/Sync/projects/NSF_AP/case_studies/"
case_dir = main_dir + "model_performance_KGI_2017_01/"
data_dir = case_dir + "models/MetUM/"
sonde_dir = case_dir + "ground_measurements/Escudero/radiosonde/"
prof_dir = case_dir + "figures/profiles/"
era_dir = case_dir + "models/ERA/"

# Files
filename = "1p5km/20170131T0000Z_Escudero_1p5km_m01s16i004_air_temperature.nc"
prof_file = "profile_20170131_18.txt"
sonde_file = "esc_sonde_dd_v2_2017013118.txt"
era_file = "era5_profile_20170131.nc"

# Output files
temp_fig = "temperature_profile.png"
temp_fig_eps = "temperature_profile.eps"

# Runtime parameters
save_figs = False

# Lat/lon of Escudero station
lat = -62.2014
lon = -58.9622
dtime = dt.datetime(2017, 1, 31, 18, 18, 0, 0)

# Load ERA5 profile and pick the closest lat/lon point
era = read_era5_profiles_interp(era_dir + era_file, lat, lon, dtime)

# Load radiosounding
sonde = pd.read_csv(sonde_dir + sonde_file, skiprows=[0, 1, 2, 4], sep="\s+")
sonde.columns = sonde.columns.str.replace(" ", "")

# Load profile
prof = pd.read_csv(prof_dir + prof_file)


# Get MetUM data

# Load air temperature from the 1.5 km run
ta = iris.load_cube(data_dir + filename, "air_temperature")

# extract data for 18Z on 31-Jan-2017
ta = ta.extract(
    iris.Constraint(time=iris.time.PartialDateTime(2017, 1, 31, 18))
)

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

# plot the MetUM temperature profile
qplt.plot(ta_at_point, ta_at_point.coord("altitude"), label="MetUM")

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


# Plot temperature from radiosonde, profile, and ERA5
plt.plot(era["temp"], era["alt"] * 1e3, label="ERA5")
plt.plot(sonde["Temp"] + 273.15, sonde["HeightMSL"], label="sonde")
plt.plot(prof["temp"], prof["alt"] * 1e3, "c.", label="profile")
plt.legend()
plt.xlabel("Temperature (K)")
plt.ylabel("Altitude (m)")
plt.axis([230, 280, 0, 6000])
plt.title("")

if save_figs:
    plt.savefig(temp_fig, format="png")
    plt.savefig(temp_fig_eps, format="eps")

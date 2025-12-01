import numpy as np
import spiceypy as spice
from datetime import datetime
from datetime import timedelta
import coord_tforms as tform

# NAIF codes
MOON = 301
EARTH = 399
JUPITER = 5
SATURN = 6
URANUS = 7
NEPTUNE = 8
SUN = 10

# center of schrodinger basin in selenographic coordinates
# https://wma.wmflabs.org/iframe.html?-75_132_0_0_en_2_en_-75_132&globe=Moon&page=Schr%C3%B6dinger%20(crater)&lang=en
SCHROD_LAT = -75 * np.pi / 180
SCHROD_LON = 132.0 * np.pi / 180

# coordinates from https://arxiv.org/abs/2508.16773
LN_LAT = -23 * np.pi/180
LN_LON = 182 * np.pi/180
class LunarObs:
    def __init__(self, observer_lat = SCHROD_LAT, observer_lon = SCHROD_LON):
        self.LAT = observer_lat
        self.LON = observer_lon
        self.OUTWARD_NORMAL = tform.lat_lon_to_rec(self.LAT, self.LON)

        self.load_kernels()

        # local coordinates here are [cos(h)cos(-A), cos(h)sin(-A), sin(h)]
        # A: azimuth measured west from south (by convention)
        # h: altitude
        # the MOON_ME coordinate axes can be rotated into the local coordinate axes via a rotation by (np.pi/2 - observer_lat) about the MOON_ME_Y axis, followed by a rotation by observer_lon about the MOON_ME_Z axis.  Hence, the matrix below transforms a vector that is expressed in the local coordinate system into the same vector expressed in the MOON_ME coordinate system
        self.local_to_moon_coords = np.matmul(tform.rotz(observer_lon), tform.roty(np.pi / 2 - observer_lat))


    @staticmethod
    def load_kernels():
        # bsp: binary spk (spacecraft and planet kernel), positions of major solar system bodies (including the moon), https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
        spice.furnsh("./spice_kernels/de440s.bsp")

        # three kernels below: moon orientation (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/)
        # bpc: binary pck (planetary constants (e.g. planet shape and orientation) kernel)
        spice.furnsh("./spice_kernels/moon_pa_de421_1900-2050.bpc")
        # tf: text frame kernel
        spice.furnsh("./spice_kernels/moon_080317.tf")
        spice.furnsh("./spice_kernels/moon_assoc_me.tf")

        # tls: text leap second kernel
        spice.furnsh("./spice_kernels/naif0012.tls")


    # get the position of 'body' relative to the moon in moon mean-earth (MOON_ME) coordinates
    @staticmethod
    def get_pos(body, times, ref = "MOON_ME"):
        body_pos = np.zeros((times.size, 3))
        for i, time in enumerate(times):
            body_pos[i] = spice.spkezp(
                targ=body, et=spice.datetime2et(time), ref=ref, abcorr="LT", obs=MOON
            )[0]
        body_pos = body_pos / ((np.sum(body_pos ** 2, axis=1)) ** (1 / 2))[:, np.newaxis]
        return body_pos


    def get_altitude(self, body, times):
        body_pos = self.get_pos(body, times)
        body_altitude = (
            90 - np.arccos(np.matmul(body_pos, self.OUTWARD_NORMAL)) * 180 / np.pi
        )
        return body_altitude

    @staticmethod
    def get_set_rise_idxs(below_horizon):
        set_idxs = below_horizon & ~np.roll(below_horizon, 1)
        rise_idxs = below_horizon & ~np.roll(below_horizon, -1)

        set_idxs[[0, -1]] = False
        rise_idxs[[0, -1]] = False

        set_idxs = np.flatnonzero(set_idxs)
        rise_idxs = np.flatnonzero(rise_idxs)

        # require full nights
        if np.min(rise_idxs) < np.min(set_idxs):
            rise_idxs = rise_idxs[1:]
        if np.max(set_idxs) > np.max(rise_idxs):
            set_idxs = set_idxs[:-1]

        return set_idxs, rise_idxs

    # print out all lunar nights that overlap with this year/month
    def print_lunar_night(self, year, month):
        utc_times = np.arange(
            datetime(year, month, 1) - timedelta(days=20),
            datetime(year, month, 1) + timedelta(days=50),
            timedelta(hours=1),
        ).astype(datetime)
        sun_altitude = self.get_altitude(SUN, utc_times)
        set_idxs, rise_idxs = self.get_set_rise_idxs(sun_altitude < 0)
        for set_time, rise_time in zip(utc_times[set_idxs], utc_times[rise_idxs]):
            if not ((set_time.month == month) or (rise_time.month) == month):
                continue
            else:
                print("night:")
                print("      start: " + str(set_time))
                print("      end: " + str(rise_time))
                print("")

    def get_lunar_nights(self, year, delta_hour_search=1, delta_min_times=15):
        utc_times = np.arange(
            datetime(year, 1, 1), datetime(year, 12, 31), timedelta(hours=1)
        ).astype(datetime)
        sun_altitude = self.get_altitude(SUN, utc_times)
        set_idxs, rise_idxs = self.get_set_rise_idxs(sun_altitude < 0)
        nights = []
        for set_time, rise_time in zip(utc_times[set_idxs], utc_times[rise_idxs]):
            nights.append(
                np.arange(set_time, rise_time, timedelta(minutes=delta_min_times)).astype(
                    datetime
                )
            )
        return nights


    # get galactic latitude and longitude of outward normal from observer coordinates
    def get_pointing(self, utc_times):
        pointing = np.zeros((utc_times.size, 3))
        for i, time in enumerate(utc_times):
            rot = spice.pxform("MOON_ME", "J2000", spice.datetime2et(time))
            pointing[i] = np.matmul(rot, self.OUTWARD_NORMAL)
        return tform.rec_to_lat_lon(*pointing.T)


    # plot the track of the outward normal on the sky over time
    def plot_day(self, utc_times, frequency, map_nside=512):
        lat, lon = self.get_pointing(utc_times)

        skymap = gsm.generate(frequency)
        skymap = hp.ud_grade(skymap, map_nside)
        skymap = np.log10(skymap)

        hp.mollview(
            skymap,
            coord="G",
            title="{} {:.0f} MHz, basemap: {}".format(gsm.name, frequency, gsm.basemap),
        )
        hp.projplot(np.pi / 2 - lat, lon, "w.", markersize = 1)






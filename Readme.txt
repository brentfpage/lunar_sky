This is a program for computing elevation angles of solar system bodies as
viewed from a given location on the lunar surface.  The heavy lifting is done
by the SPICE toolkit from NASA/JPL/NAIF
(https://naif.jpl.nasa.gov/naif/toolkit.html), wrapped by SpiceyPy
(https://pypi.org/project/spiceypy/). An observer object is created via, e.g.,
lunar_observer = LunarObs(LN_LAT, LN_LON) where LN_LAT and LN_LON are an
observer's latitude and longitude in selenographic coordinates.  Functions such
as lunar_observer.get_altitude(...) can then be called that are specific to
this observer's latitude and longitude.  This function get_altitude(...)
returns the angle between the horizon and the line of sight to some solar
system body.  There are other functions, such as get_pos(...), that involve
moon-related ephemeris generally and aren't specific to a certain viewing
location.  These are still called using an observer object, as in
lunar_observer.get_pos(...) but the result does not depend on the properties of
lunar_observer.

The program's basic functionality is demonstrated in driver.py, which can be
called from the command line as python3 driver.py

Prerequisities: most notably, spiceypy, which can be downloaded via pip3
install spiceypy --break-system-packages where the --break-system-packages flag
should be added at the user's own risk but in my experience is required and
harmless.

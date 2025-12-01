import lunar_sky as ls
import datetime
import numpy as np

schrod_obs = ls.LunarObs() # observer from schrodinger basin
saturn_alt_from_schrod = schrod_obs.get_altitude(ls.SATURN,np.array([datetime.datetime(2025,9,20)]))

ln_obs = ls.LunarObs(ls.LN_LAT,ls.LN_LON) # latitude and longitude listed in https://arxiv.org/abs/2508.16773
saturn_alt_from_ln = ln_obs.get_altitude(ls.SATURN,np.array([datetime.datetime(2025,9,20)]))

print(saturn_alt_from_schrod)
print(saturn_alt_from_ln)
print(schrod_obs.get_pos(ls.SATURN,np.array([datetime.datetime(2025,9,20)])))
print(ln_obs.get_pos(ls.SATURN,np.array([datetime.datetime(2025,9,20)])))

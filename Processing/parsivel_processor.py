import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 

import glob
import pydsd
import xarray as xr
import os
from pydsd.io import ParsivelReader as pr
import datetime
import numpy as np


oldfiles = sorted(glob.glob('/data/accp/a/snesbitt/scamp/parsivel/*.MIS'))

for ifile in oldfiles:

    dsd = pydsd.read_parsivel(ifile)
    
    time_val = [datetime.datetime.fromtimestamp(x) + datetime.timedelta(hours=12) for x in dsd.time["data"]]

    yyyymmdd = time_val[0].strftime('%Y%m%d')
    yyyymmddhhmmss = time_val[0].strftime('%Y%m%d_%H%M%S')
    
    if os.path.isfile('/data/accp/a/snesbitt/scamp/parsivel/nc/'+yyyymmdd+'/'+yyyymmddhhmmss+'.nc') != True:
        print('processing '+yyyymmddhhmmss)
    
        raw = pr.ParsivelReader(ifile)
#        velocity = np.array(raw.velocity["data"],dtype='float64')

        data = xr.Dataset(data_vars=dict(
             Nd=(["time","diameter"], dsd.fields["Nd"]["data"]),
             num_particles=(["time"], dsd.fields['num_particles']["data"]),
             dBZ=(["time"], raw.Z),
             rr=(["time"], raw.rain_rate),
             spectrum=(["time", "velocity","diameter"], np.reshape(raw.raw, (1,32,32))),
         ),
         coords=dict(
             time=(["time"], time_val),
             diameter=(["diameter"], dsd.diameter["data"]),
             velocity=(["velocity"], raw.velocity["data"]),
             bin_edges=(["bin_edges"], dsd.bin_edges["data"]),
             v_spread=(["v_spread"], raw.v_spread)
         ),
         attrs=dict(description="SCAMP Parsivel Data"))
        yyyymmdd = time_val[0].strftime('%Y%m%d')
        yyyymmddhhmmss = time_val[0].strftime('%Y%m%d_%H%M%S')
        os.system('mkdir -p /data/accp/a/snesbitt/scamp/parsivel/nc/'+yyyymmdd)
        data.to_netcdf('/data/accp/a/snesbitt/scamp/parsivel/nc/'+yyyymmdd+'/'+yyyymmddhhmmss+'.nc')



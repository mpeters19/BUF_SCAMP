import numpy as np
import xarray as xr
from netCDF4 import Dataset
import warnings
from datetime import date
from tqdm import tqdm
from scipy.optimize import curve_fit
import os

def calc_dsd(date):

    warnings.filterwarnings("ignore", category=RuntimeWarning)
    ds = xr.open_mfdataset('/data/accp/a/mp46/Parsivel_data/DSD/'+date+'_DSD_BUF.nc')

    gauge_and_wx_station_data = xr.open_mfdataset('/data/accp/a/mp46/gauge_wx_sensor_filtered_data_BUF.nc')

    dsd = ParsivelPSDBUF(ds)

    dsd.get_precip_params(gauge_and_wx_station_data)

    dsd.saveDSD_as_nc(ds,date)

    ds.close()
    gauge_and_wx_station_data.close()
    #os.remove('/data/accp/a/mp46/Parsivel_data/'+date+'_ND_BUF.nc')
    return dsd

def power_law(x, a, b):
        return a*np.power(x, b)

class ParsivelPSDBUF(object):

    def __init__(self,ds):
        #some data stored as tuples if applicable
        #first element = average or total across all 32 DSD bins
        #second element = value for each bin individually (plus time-dimension)
        self.ds = ds #initializing ds as the original disdrometer data file
        self.time = self.ds.time 
        self.timedim = len(ds.time) #time dimension
        #DSD is the number of drops per volume of air for a given drop size interval (i.e. for the 32 bins)
        self.terminal_velocities = np.zeros((len(ds.time),32))
        self.Vt_parameters = np.zeros((len(ds.time),2))
        self.m_D_coefficient = np.zeros(self.timedim)

        print('Variables Initialized')

    def get_precip_params(self, gauge_and_wx_station_data):
        
        print('Calculations in progess')
        
        time_ds1 = self.time.values
        time_ds2 = gauge_and_wx_station_data.UTC.values

        # Find the common timestamps
        common_timestamps = np.intersect1d(time_ds1, time_ds2)

        # Find the indices of common timestamps in ds2
        common_indices_ds2 = np.where(np.isin(time_ds2, common_timestamps))[0]
        
        print('Common time indicies found')
        
        for td in list(range(0,len(self.ds.time))):
            print(td)
           
            snow_dia_vel = []

            for i in range(0, len(self.drop_diameter)):
                local_dia = self.drop_diameter[i]
                for j in range(0, len(self.v_parsivel)):
                    local_veloc = self.v_parsivel[j]
                    count = self.ds.spectrum[td, i, j]
                    try:
                        for _ in range(int(count)):
                            snow_dia_vel.append((local_dia, local_veloc))
                    except ValueError:
                        snow_dia_vel = []
                        
            snow_triple_time_dia_vel = np.array(snow_dia_vel)

            if len(snow_triple_time_dia_vel) == 0:
                current_vt = np.zeros(32)
                self.Vt_parameters[td,0] = 0
                self.Vt_parameters[td,1] = 0
            else:
                snow_pars, snow_cov = curve_fit(f=power_law, xdata=snow_triple_time_dia_vel[:,0]/1000, ydata=snow_triple_time_dia_vel[:,1], p0=[0.006, 0.05], bounds=(0, [np.inf, 1]), check_finite=False)
                current_vt = (snow_pars[0])*np.power(np.array(self.drop_diameter)/1000, snow_pars[1])
                self.Vt_parameters[td,0] = snow_pars[0]
                self.Vt_parameters[td,1] = snow_pars[1]

            self.terminal_velocities[td,:] = current_vt
        
        print('Calulating a')
        gauge_preciprate_oneminute = gauge_and_wx_station_data.Processed_precip_rate_mmhr[common_indices_ds2[0]:(common_indices_ds2[-1]+1)]

        rolling = gauge_preciprate_oneminute.rolling(UTC=30, center=True)
        gauge_data_rolling = rolling.mean()

        gauge_preciprate = np.array(np.divide(gauge_data_rolling, 3.6)) # to convert from mm/hr = kg m^-2 hr^-1 to g m^-2 s^-1
        diam_squared = np.power(np.divide(np.array(self.drop_diameter),10),2)
        denom = np.array(np.sum(self.ds.dsd * self.drop_spread * self.terminal_velocities * diam_squared, axis = 1))

        self.m_D_coefficient = np.divide(gauge_preciprate, denom)
         
        print('Calculations done')

    def saveDSD_as_nc(self,ds,date):
        
        print('Saving netCDF file for '+date+'...')
        ds["terminal_velocities"] = (["time", "diameter"], self.terminal_velocities.data)
        #ds["Vt_coefficient"] = (["time"], self.Vt_parameters[0].data)
        #ds["Vt_exponent"] = (["time"], self.Vt_parameters[1].data)
        ds["m_D_coefficient"] = (["time"], self.m_D_coefficient.data)
        
        ds.to_netcdf('/data/accp/a/mp46/Parsivel_data/' + date + '_Vt_mD_BUF.nc')  # saving netCDF to specified path
        print("NetCDF file saved successfully.")

    drop_diameter = [
        0.062, 0.187, 0.312, 0.437, 0.562, 0.687, 0.812, 0.937, 1.062, 1.187, 1.375, 1.625,
        1.875, 2.125, 2.375, 2.750, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 
        13, 15, 17, 19, 21.5, 24.5] #diameters from the OTT Parsivel2 manual in mm

    drop_spread = [
        0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.250,
        0.250, 0.250, 0.250, 0.250, 0.500, 0.500, 0.500, 0.500, 0.500, 1.000, 1.000,
        1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 2.000, 2.000, 3.000, 3.000] #also delta

    v_parsivel = [
        0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.1, 1.3, 1.5, 1.7, 1.9,
        2.2, 2.6, 3, 3.4, 3.8, 4.4, 5.2, 6.0, 6.8, 7.6, 8.8, 10.4, 12.0, 13.6, 15.2,
        17.6, 20.8]

    
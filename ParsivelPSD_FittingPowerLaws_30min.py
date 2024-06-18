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

            # Example data
            # velocity: 32, diameter: 32, data: 2D matrix of counts
            velocity_bins = self.v_parsivel  # Example velocity bins
            diameter_bins = self.drop_diameter  # Example diameter bins
            counts_matrix = self.ds.spectrum[td].values # Example counts matrix

            # Prepare the data for fitting
            velocities, diameters = np.meshgrid(velocity_bins, diameter_bins, indexing='ij')
            velocities_flat = velocities.flatten()
            diameters_flat = diameters.flatten()
            counts_flat = counts_matrix.flatten()

            # Filter out zero counts
            mask = counts_flat > 0
            filtered_velocities = velocities_flat[mask]
            filtered_diameters = diameters_flat[mask]
            filtered_counts = counts_flat[mask]

            # Prepare the data for curve fitting
            xdata = filtered_diameters / 1000  # Convert to meters if needed
            ydata = filtered_velocities

            if ydata.size == 0:
                a, b = 0, 0
            else:
                # Fit the power law
                popt, pcov = curve_fit(power_law, xdata, ydata, p0=[0.006, 0.05], bounds=(0, [np.inf, 1]), check_finite=False)
                # Extract fitted parameters
                a, b = popt

            # Calculate fitted velocities for the original diameter bins
            fitted_velocities = power_law(np.array(diameter_bins) / 1000, a, b)
            self.terminal_velocities[td,:] = fitted_velocities
           
            self.Vt_parameters[td,0] = a
            self.Vt_parameters[td,1] = b

            # Calculate weighted sum and mean terminal velocity
            spectrum_sum = np.sum(self.ds.spectrum[td, :, :])
            weighted_sum = np.sum(fitted_velocities * np.sum(self.ds.spectrum[td, :, :], axis=0), axis=0)
            mean_terminal_velocity = np.where(spectrum_sum != 0, weighted_sum / spectrum_sum, 0)
            self.mean_fitted_vt[td] = mean_terminal_velocity

            spectrum_sum = np.sum(np.sum(self.ds.spectrum[td, :, :], axis=1))
            if spectrum_sum != 0:
                self.mean_measured_v[td] = np.sum(np.sum(self.ds.spectrum[td, :, :], axis=1) * self.v_parsivel) / spectrum_sum
            else:
                self.mean_measured_v[td] = 0
        
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
        ds["terminal_velocities"] = (["time", "diameter"], self.terminal_velocities)
        ds["Vt_coefficient"] = (["time"], self.Vt_parameters[:, 0])
        ds["Vt_exponent"] = (["time"], self.Vt_parameters[:, 1])
        ds["m_D_coefficient"] = (["time"], self.m_D_coefficient)
        ds["mean_fitted_vt"] = (["time"], self.mean_fitted_vt)
        ds["mean_measured_v"] = (["time"], self.mean_measured_v)
        
        ds.to_netcdf('/data/accp/a/mp46/Parsivel_data/' + date + '_Vt_mD_BUF30min_overhaul.nc')  # saving netCDF to specified path
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

    

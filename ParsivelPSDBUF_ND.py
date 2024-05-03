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

    # Takes the date you want processed as an input, and then opens the files from that day (since they are hourly files)

    ds_10s = xr.open_mfdataset('/data/accp/a/snesbitt/scamp/parsivel/nc_daily/'+date+'_*_SCAMP_Parsivel.nc')
    ds = ds_10s.resample(time='1MIN').sum() #data is originally recorded every 10s, this resamples it to 1 min by summing all data within a minute

    dsd = ParsivelPSDBUF(ds)

    dsd.get_precip_params()

    dsd.saveDSD_as_nc(date)

    ds_10s.close()
    ds.close()
    return dsd

class ParsivelPSDBUF(object):

    '''
    ParsivelPSDBUF takes a ProcessParsivel instance and computes the DSD using the methods from Tokay et al. (2014). 
    This code is based on the get_rainparams.pro code written by Ali Tokay and converted to IDL by Dave Wolff. 

    The reason for keeping the ParsivelDSD as a separate class is that the user may be interested in seeing 
    characteristics of the DSD that do not appear in the final processed dataset...for instance, the contribution
    to LWC from each Parsivel bin. By using this class, one can branch off from various steps on the processing 
    process and make plots or analyze the data.

    In some cases, the data is saved as a tuple, with the contribution to rainrate from each drop bin range 
    included as the second element in the tuple. 

    Notes:
    Currently uses the corrected fall velocities from Ali Tokay...changing the assumed fall velocities will 
    affect the drop parameters. 
    '''

    def __init__(self,ds):
        #some data stored as tuples if applicable
        #first element = average or total across all 32 DSD bins
        #second element = value for each bin individually (plus time-dimension)
        self.ds = ds #initializing ds as the original disdrometer data file
        self.time = self.ds.time 
        self.timedim = len(ds.time) #time dimension
        #DSD is the number of drops per volume of air for a given drop size interval (i.e. for the 32 bins)
        self.dsd = np.zeros((self.timedim,32)) #dsd 

    def get_precip_params(self):
        #Takes a processed_parsivel instance and calculates N(D)
        delta_t = 60 # units: seconds
        p2_area = 180.*(30.-(np.array(self.drop_diameter)/2.)) # same as effective_sampling_area in PyDSD, parsivel laser volume is 180x30x1 mm, (subtracting off partial drops on edge)

        for td in list(range(0,len(self.ds.time))):
            
            #N(D), Tokay et al. 2014, Eq 6, note: denom is in per m^3*mmbin, #10^6 converts to m^3*mm instead of mm^3*m
            self.dsd[td,:] = 1e6 * np.dot(np.swapaxes(self.ds.spectrum[td], 0, 1), 1 / np.array(self.v_parsivel),)/ (p2_area * np.array(self.drop_spread) * delta_t)
            
    def saveDSD_as_nc(self,date):
        
        '''Takes the dsd object and saves it as a netCDF file'''
        
        data = xr.Dataset(data_vars=dict(
             dsd = (["time","diameter"], self.dsd.data), #The number of drops per volume of air for a given drop size interval (ie for the 32 bins)
             spectrum = (["time", "velocity","diameter"], self.ds.spectrum.data),
             
         ),
         coords=dict(
             time=(["time"], self.time.data), #assigning time as a coordinate
             diameter=(["diameter"], self.drop_diameter), #assigning drop diameter bins as a coordinate
             velocity=(["velocity"], self.v_parsivel), #assigning velocity bins as a coordinate
         ),
         attrs=dict(description="BUF Parsivel Data")) # adding a description for netCDF file

        print('Saving netCDF file for '+date+'...')
        data.to_netcdf('/data/accp/a/mp46/Parsivel_data/' + date + '_DSD_BUF.nc')  # saving netCDF to specified path
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
    

    
import numpy as np
import pdb
import datetime
import scipy.stats.mstats
import glob
import zipfile
from math import gamma
import math
import xarray as xr
from netCDF4 import Dataset
from utility import filter
import pytmatrix
import scipy
from scipy.optimize import curve_fit
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, radar, tmatrix_aux, refractive
from datetime import date
import time
from tqdm import tqdm

def calc_dsd(date):
    '''
    Funtion that wraps all three Parsivel classes together to make the
    final DSD object

    Use to get the data to make plots
    '''
    ds_10s = xr.open_mfdataset('/data/accp/a/snesbitt/scamp/parsivel/nc_daily/'+date+'_*_SCAMP_Parsivel.nc')
    ds = ds_10s.resample(time='1MIN').mean()

    dsd = ParsivelPSDBUF(ds)

    dsd.get_precip_params()

    dsd.saveDSD_as_nc(date)
    return dsd

class ParsivelPSDBUF(object):

    '''
    ParsivelDSD takes a ProcessParsivel instance and computes the DSD using
    the methods from Tokay et al. (2014). This code is based on the
    get_rainparams.pro code written by Ali Tokay and converted to IDL by 
    Dave Wolff. 

    The reason for keeping the ParsivelDSD as a separate class is that the
    user may be interested in seeing characteristics of the DSD that do
    not appear in the final processed dataset...for instance, the contribution
    to LWC from each Parsivel bin. By using this class, one can branch off from
    various steps on the processing process and make plots or analyze the data.

    In some cases, the data is saved as a tuple, with the contribution to rainrate
    from each drop bin range included as the second element in the tuple. 

    The PlotParsivel class can be used to make plots of this data.

    Notes:
    Currently uses the corrected fall velocities from Ali Tokay...changing the 
    assumed fall velocities will affect the drop parameters. 

    Future additions:
    Standard deviation of drop size
    '''

    def __init__(self,ds):
        #some data stored as tuples if applicable
        #first element = average or total across all 32 DSD bins
        #second element = value for each bin individually (plus time-dimension)
        self.ds = ds
        self.time = self.ds.time 
        timedim = len(ds.time) #time dimension
        #ndrops is the number of drops (non-normalzed) in each volume of air
        #DSD is the number of drops per volume of air for a given drop size interval (i.e. for the 32 bins)
        self.dsd = np.zeros((timedim,32)) #dsd 
        self.vvd = np.zeros((timedim,32))
        self.mass_warm = (np.zeros(timedim),np.zeros((timedim,32)))
        self.mass_cold = (np.zeros(timedim),np.zeros((timedim,32)))
        self.mass_conv = (np.zeros(timedim),np.zeros((timedim,32)))
        self.IWC = (np.zeros(timedim),np.zeros((timedim,32)))
        self.LWC = (np.zeros(timedim),np.zeros((timedim,32)))
        self.liquid_mass = (np.zeros(timedim),np.zeros((timedim,32)))
        #self.parsivel_velocity
        #drop_conc is the number of drops per volume of air. It is NOT normalized by drop size interval like DSD
        self.z = (np.zeros(timedim),np.zeros((timedim,32))) #reflectivity factor from each drop size bin
        self.ze = (np.zeros(timedim),np.zeros((timedim,32))) #reflectivity factor from each drop size bin
        self.dbz = np.zeros(timedim)
        self.rainrate = (np.zeros(timedim),np.zeros((timedim,32,32))) #rainrate from each drop size bin
        self.dm_moments = np.zeros(timedim) #mass-weighted mean diameter
        self.moments = np.zeros((timedim,8)) #drop moments


    def get_precip_params(self):
        
        #Takes a processed_parsivel instance and calculates the precip params
         
        # parsivel_matrix: 1ength 1024 matrix of Parsivel diam/fspd
        # timerain: time of period in seconds #self.proc_p2.time_interval
        # wxcode: needed to differentiate between rain and snow
        # Note: if frozen precipitation detected, no rain rate or LWC is returned
        
        #add loop here to go through each time dimension
        #timerain has to be calculated individually for each record in case data is missing 
        #what usually happens is that maybe 1 min of data per day is missing in 10-30s intervals
        timerain = np.array(1) #minute
        for td in list(range(0,len(self.ds.time))):
            #get correct time multiplier
            time_mult = 60  #units: seconds
            #time_div = 60 #units: s/min
            matrix = np.array(self.ds.spectrum[td,:,:])
            
            #ndrops = np.sum(matrix)
            #moments: 0) concen, 1) mean diam, 2) surface area conc, 3)lwc, 6)z
            #self.moments[td,3] = 0.0 #for moments
            #self.moments[td,4] = 0.0 #for moments
            #self.moments[td,2] = 0.0 #for moments
            #self.moments[td,6] = 0.0 #for moments
            #Process each drop bin (diameter, velocity):
            for dind,dbin in enumerate(self.drop_diameter):
                for vind,vbin in enumerate(self.v_parsivel):
                    drops = matrix[vind,dind]
                    #p2_area = filter.parsivel_sampling_area(dind)
                    p2_area = 180.*(30.-(dbin/2.))/100.
                    #denominators
                    denominator = time_mult * p2_area * vbin 
                    #denom2 = time_mult * p2_area * vbin * self.drop_spread[dind] #per m^3*mmbin
                    denom2 = time_mult * p2_area * vbin * self.drop_spread[dind] * 100 #per m^3*mmbin

                    self.dsd[td,dind] += (1.e6 * drops)/denom2 #10^6 converts to m^3*mm instead of mm^3*m

                    moment3 = np.power(self.drop_diameter[dind], 3)
                    self.moments[td,3] += np.multiply(np.multiply(moment3, self.dsd[td,dind]), self.drop_spread[dind])
                    moment4 = np.power(self.drop_diameter[dind], 4)
                    self.moments[td,4] += np.multiply(np.multiply(moment4, self.dsd[td,dind]), self.drop_spread[dind])
                    moment6 = np.power(self.drop_diameter[dind], 6)
                    self.moments[td,6] += np.multiply(np.multiply(moment6, self.dsd[td,dind]), self.drop_spread[dind])

                    # #reflectivity factor (6th power of diameter, normalized for area, time)
                    self.z[1][td,dind] += drops * 1.e6 * dbin**6/denominator
                    # equivalent reflectivity, Smith 1984
                    self.ze[1][td,dind] += (0.176/0.93)*self.z[1][td,dind]
                    #reflectivity weighted vvd
                    self.vvd[td,vind] = np.divide(moment6*vbin, moment6)
                    # using PyDSD
                    #self.rainrate[1][td,dind] += 0.6 * np.pi * 1e-03 * np.dot(np.multiply(self.proc_p2.rainrate[td],np.multiply(self.ndrops[0][td],self.drop_spread[dind])),np.array(self.drop_diameter[dind]) ** 3)

            #Williams et al. 2014, mass-weighted mean diameter = ratio of 4th to 3rd moments applicable for rain 
            self.dm_moments[td] = self.moments[td,4] / self.moments[td,3] 

            # ice mass, Heymsfield 2010
            self.mass_spectrum_warm[1][td] = np.multiply(self.dsd[td],self.heymsfield_warm_mass)
            self.mass_spectrum_cold[1][td] = np.multiply(self.dsd[td],self.heymsfield_cold_mass)       
            self.mass_spectrum_conv[1][td] = np.multiply(self.dsd[td],self.heymsfield_conv_mass)       
            self.mass_spectrum_warm[0][td] = np.sum(self.mass_spectrum_warm[1][td,:])
            self.mass_spectrum_cold[0][td] = np.sum(self.mass_spectrum_cold[1][td,:])
            self.mass_spectrum_conv[0][td] = np.sum(self.mass_spectrum_conv[1][td,:])

            # IWC, Heymsfield 2005
            #mass_multiplication = np.multiply(np.divide(self.dsd[td],10^5),self.heymsfield_warm_mass)
            self.IWC[1][td] = np.multiply(10^6,self.dsd[td],self.heymsfield_warm_mass,self.drop_spread)
            #self.IWC[1][td] = np.multiply(self.mass_warm[1][td],self.drop_spread)
            self.IWC[0][td] = np.sum(self.IWC[1][td,:])

            # LWC, lecture 1 slide 41 of 510
            moments_for_LWC = np.power(np.divide(self.drop_diameter,1000),3)
            coefficient = np.divide(np.pi*1000000,6) #density here is 1000 kg/m^3 multiplied by a conversion factor to get kg->g
            following_sigma = np.multiply(moments_for_LWC,self.dsd[td])
            self.LWC[1][td] = np.multiply(coefficient,following_sigma)
            self.LWC[0][td] = np.sum(self.LWC[1][td,:])

            # liquid mass, Williams et al. 2014
            moments_for_liquid_mass = np.power(np.divide(self.drop_diameter,10),3)
            coefficient = np.divide(np.pi*1,6*1000) #density here is 1 g/cm^3 
            diameter_calc = np.multiply(moments_for_liquid_mass,self.dsd[td])
            self.liquid_mass[1][td] = np.multiply(coefficient,diameter_calc)
            self.liquid_mass[0][td] = np.sum(self.LWC[1][td,:])

            # equivalent reflectivity, Smith 1984
            self.ze[0][td] = np.sum(self.ze[1][td,:])
            # #reflectivity factor (6th power of diameter, normalized for area, time)
            self.z[0][td] = np.sum(self.z[1][td,:])
            if self.z[0][td] > 0:
                self.dbz[td] = 10 * np.log10(self.z[0][td]) #z to dbz
            else:
                self.dbz[td] = float('nan')
        '''
        summed_spectrum = self.ds.spectrum.sum(axis=0)
        separate_diameters = []
        separate_velocities = []

        for i in tqdm(range(len(summed_spectrum))):
            for j in range(len(summed_spectrum)):
                #niftycount = math.floor(summed_spectrum[j,i].values)
                niftycount = summed_spectrum[j,i].values//1
                if niftycount > 0:
                    separate_diameters.extend([float(self.drop_diameter[i])]*niftycount)
                    separate_velocities.extend([float(self.v_parsivel[j])]*niftycount)

        # Fit the dummy power-law data
        self.pars, cov = curve_fit(f=power_law, xdata=separate_diameters, ydata=separate_velocities, p0=[0, 0], bounds=(-np.inf, np.inf))
        # Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
        #stdevs = np.sqrt(np.diag(cov))
        # Calculate the residuals
        #res = ds.velocity - power_law(ds.diameter, *pars)
        def power_law(x, a, b):
            return a*np.power(x, b)
        '''
    def saveDSD_as_nc(self,date):

        '''Takes the dsd object and saves it as a netCDF file'''
        data = xr.Dataset(data_vars=dict(
             #drop_concentration_bins = (["time","diameter"], self.drop_conc[1]), #The number of drops per volume of air, value for each bin individually (plus time-dimension)
             #total_drop_concentration = (["time"],self.drop_conc[0]), #The number of drops per volume of air, average or total across all 32 DSD bins
             avg_rainrate = (["time"], self.ds.rr), #mean Parsivel rain rate value recorded in each data chunk
             avg_dbz = (["time"], self.ds.dBZ), #mean Parsivel dBZ recorded in each data chunk
             #temperature = (["time"], self.proc_p2.temperature), #The mean parsivel temperature value of each data chunk
             mean_diameter = (["time"], self.dm_moments[:]), #Mass-weighted mean diameter
             #max_diameter = (["time"], self.dmax), #Max drop size
             #sigma_mass = (["time"], self.sigma_m), #Variance of mass spectrum
             calc_dBZ = (["time"], self.dbz[:]), #mean Parsivel rain rate value recorded in each data chunk
             #number_drops_bins = (["time","diameter"], self.ndrops[1]), #Number of drops non-normalized in each volume of air, value for each bin individually (plus time-dimension)
             #total_number_drops = (["time"],self.ndrops[0]), #Number of drops non-normalized in each volume of air, average or total across all 32 DSD bins
             #rain_accumulation_bins = (["time","diameter"], self.rainaccum[1]), #rain accumulation from each bin, value for each bin individually (plus time-dimension)
             #total_rain_accumulation = (["time"], self.rainaccum[0]), #rain accumulation, average or total across all 32 DSD bins
             original_dsd = (["time","diameter"], self.ds.Nd),
             dsd = (["time","diameter"], self.dsd), #The number of drops per volume of air for a given drop size interval (ie for the 32 bins)
             vvd = (["time","velocity"], self.vvd),
             Heymsfield_mass_warm_bins = (["time","diameter"], self.mass_warm[1]),
             Heymsfield_mass_warm_total = (["time"], self.mass_warm[0]),
             Heymsfield_mass_warm_no_dsd = (["diameter"], self.heymsfield_warm_mass),
             Heymsfield_mass_cold_bins = (["time","diameter"], self.mass_cold[1]),
             Heymsfield_mass_cold_total = (["time"], self.mass_cold[0]),
             Heymsfield_mass_conv_bins = (["time","diameter"], self.mass_conv[1]),
             Heymsfield_mass_conv_total = (["time"], self.mass_conv[0]),
             IWC_bins = (["time","diameter"], self.IWC[1]),
             IWC_total = (["time"], self.IWC[0]),
             LWC_bins = (["time","diameter"], self.LWC[1]),
             LWC_total = (["time"], self.LWC[0]),
             raindrop_mass_bins = (["time","diameter"], self.liquid_mass[1]),
             raindrop_mass_total = (["time"], self.liquid_mass[0]),
             z_bins = (["time","diameter"], self.z[1]), #Reflectivity factor from each drop size bin, value for each bin individually (plus time-dimension)
             total_z = (["time"], self.z[0]), # reflectivity factor, average or total across all 32 DSD bins
             ze_bins = (["time","diameter"], self.ze[1]), #Reflectivity factor from each drop size bin, value for each bin individually (plus time-dimension)
             total_ze = (["time"], self.ze[0]), # reflectivity factor, average or total across all 32 DSD bins
             #rain_rate_bins = (["time","velocity","diameter"], self.rainrate[1]), #rain accumulation from each bin, value for each bin individually (plus time-dimension)
             #total_rainrate = (["time"], self.rainrate[0]), #rain accumulation, average or total across all 32 DSD bins
             
             #spectrum = (["time", "velocity","diameter"], np.reshape(self.ds.spectrum, (1440,32,32)))
             spectrum = (["time", "velocity","diameter"], self.ds.spectrum)
         ),
         coords=dict(
             time=(["time"], self.time), #assigning time as a coordinate
             diameter=(["diameter"], self.drop_diameter), #assigning drop diameter bins as a coordinate
             velocity=(["velocity"], self.v_parsivel), #assigning velocity bins as a coordinate
         ),
         attrs=dict(description="BUF Parsivel Data")) # adding a description for netCDF file
        #yyyymmdd = self.time.dt.strftime('%Y%m%d')
        #yyyymmdd = str(yyyymmdd)
        data.to_netcdf('/data/keeling/a/mp46/Research/Processed_Data/'+date+'_BUFtest.nc') #saving netCDF to specified path

    
    drop_diameter = [
        0.062, 0.187, 0.312, 0.437, 0.562, 0.687, 0.812, 0.937, 1.062, 1.187, 1.375, 1.625,
        1.875, 2.125, 2.375, 2.750, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 
        13, 15, 17, 19, 21.5, 24.5] #diameters from the OTT Parsivel2 manual

    drop_spread = [
        0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.250,
        0.250, 0.250, 0.250, 0.250, 0.500, 0.500, 0.500, 0.500, 0.500, 1.000, 1.000,
        1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 2.000, 2.000, 3.000, 3.000] #also delta

    v_parsivel = [
        0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.1, 1.3, 1.5, 1.7, 1.9,
        2.2, 2.6, 3, 3.4, 3.8, 4.4, 5.2, 6.0, 6.8, 7.6, 8.8, 10.4, 12.0, 13.6, 15.2,
        17.6, 20.8]
    
    heymsfield_power = np.power(np.divide(drop_diameter,10),2.1)
    heymsfield_warm_mass = np.multiply(0.00359,heymsfield_power)
    heymsfield_cold_mass = np.multiply(0.00574,heymsfield_power)
    heymsfield_conv_mass = np.multiply(0.00630,heymsfield_power)



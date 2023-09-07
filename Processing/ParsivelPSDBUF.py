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
from utility import configuration
import pytmatrix
import scipy
from scipy.optimize import curve_fit
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, radar, tmatrix_aux, refractive
from datetime import date
import time
from tqdm import tqdm
from PyDSD.pydsd.DropSizeDistribution import DropSizeDistribution
from PyDSD.pydsd.io import common

def calc_dsd(date):
    '''
    Funtion that wraps all three Parsivel classes together to make the
    final DSD object

    Use to get the data to make plots
    '''
    ds_10s = xr.open_mfdataset('/data/accp/a/snesbitt/scamp/parsivel/nc_daily/'+date+'_*_SCAMP_Parsivel.nc')
    ds = ds_10s.resample(time='1MIN').sum()

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
        self.timedim = len(ds.time) #time dimension
        #ndrops is the number of drops (non-normalzed) in each volume of air
        #DSD is the number of drops per volume of air for a given drop size interval (i.e. for the 32 bins)
        self.dsd = np.zeros((self.timedim,32)) #dsd 
        self.vvd = np.zeros(self.timedim)
        self.mass_spectrum_warm = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.mass_spectrum_cold = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.mass_spectrum_conv = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.IWC = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.LWC = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.lwc = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.Williams_LWC = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.liquid_mass = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.liquid_mass_spectrum = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.liquid_mass_simple = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.liquid_mass_spectrum_simple = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.mass_spectrum_matrosov = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.z = (np.zeros(self.timedim),np.zeros((self.timedim,32))) #reflectivity factor from each drop size bin
        self.ze = (np.zeros(self.timedim),np.zeros((self.timedim,32))) #reflectivity factor from each drop size bin
        self.dbz = np.zeros(self.timedim)
        self.rainrate = (np.zeros(self.timedim),np.zeros((self.timedim,32,32))) #rainrate from each drop size bin
        self.dm_moments = np.zeros(self.timedim) #mass-weighted mean diameter
        self.moments = np.zeros((self.timedim,8)) #drop moments 
        self.D0 = np.zeros(self.timedim)
        #self.Dmax = np.zeros(self.timedim)
        self.Dm = np.zeros(self.timedim)
        self.Nt = np.zeros(self.timedim)
        self.Nw = np.zeros(self.timedim)
        self.N0 = np.zeros(self.timedim)
        self.mu = np.zeros(1)
        self.Lambda = np.zeros(1)

    def get_precip_params(self):
        
        #Takes a processed_parsivel instance and calculates the precip params
         
        # parsivel_matrix: 1ength 1024 matrix of Parsivel diam/fspd
        # timerain: time of period in seconds #self.proc_p2.time_interval
        # wxcode: needed to differentiate between rain and snow
        # Note: if frozen precipitation detected, no rain rate or LWC is returned
        
        #add loop here to go through each time dimension
        #timerain has to be calculated individually for each record in case data is missing 
        #what usually happens is that maybe 1 min of data per day is missing in 10-30s intervals
        for td in list(range(0,len(self.ds.time))):
            #get correct time multiplier
            time_mult = 60  #units: seconds
            #time_div = 60 #units: s/min
            matrix = np.array(self.ds.spectrum[td,:,:])
            vvd_num = 0
            vvd_denom = 0
            
            #ndrops = np.sum(matrix)
            #Process each drop bin (diameter, velocity):
            for dind,dbin in enumerate(self.drop_diameter):
                for vind,vbin in enumerate(self.v_parsivel):
                    drops = matrix[vind,dind]
                    #Next step uses equation (6) from Tokay et al. (2014)
                    #parsivel laser area, units: mm^2 (subtracting off partial drops on edge)

                    p2_area = 180.*(30.-(dbin/2.))
                    #denominators
                    denominator = time_mult * p2_area * vbin 
                    #denom2 = time_mult * p2_area * vbin * self.drop_spread[dind] #per m^3*mmbin
                    denom2 = time_mult * p2_area * vbin * self.drop_spread[dind] #per m^3*mmbin

                    self.dsd[td,dind] += (1.e6 * drops)/denom2 #10^6 converts to m^3*mm instead of mm^3*m

                    moment3 = np.power(self.drop_diameter[dind], 3)
                    self.moments[td,3] += np.multiply(np.multiply(moment3, self.dsd[td,dind]), self.drop_spread[dind])
                    moment4 = np.power(self.drop_diameter[dind], 4)
                    self.moments[td,4] += np.multiply(np.multiply(moment4, self.dsd[td,dind]), self.drop_spread[dind])
                    moment6 = np.power(self.drop_diameter[dind], 6)
                    self.moments[td,6] += np.multiply(np.multiply(moment6, self.dsd[td,dind]), self.drop_spread[dind])
                    
                    vol = np.pi*dbin**3/6 #volume of 1 drop in this size bin
                    self.lwc[1][td,dind] += drops*vol*1.e3/denominator #units: g/m^3 per mm bin (rho=1000 g/m^3)
                    #reflectivity factor (6th power of diameter, normalized for area, time)
                    self.z[1][td,dind] += drops * 1.e6 * dbin**6/denominator
                    # equivalent reflectivity, Smith 1984
                    self.ze[1][td,dind] += (0.176/0.93)*self.z[1][td,dind]
                    #reflectivity weighted vvd
                    vvd_num += np.sum(moment6 * (1.e6 * drops)/denom2 * self.drop_spread[dind] * vbin)
                    vvd_denom += np.sum(moment6 * (1.e6 * drops)/denom2 * self.drop_spread[dind])
                    # using PyDSD
                    #self.rainrate[1][td,dind] += 0.6 * np.pi * 1e-03 * np.dot(np.multiply(self.proc_p2.rainrate[td],np.multiply(self.ndrops[0][td],self.drop_spread[dind])),np.array(self.drop_diameter[dind]) ** 3)

            #Williams et al. 2014, mass-weighted mean diameter = ratio of 4th to 3rd moments applicable for rain 
            self.dm_moments[td] = self.moments[td,4] / self.moments[td,3] 

            # ice mass, Heymsfield 2010
            self.mass_spectrum_warm[1][td] = np.multiply(self.dsd[td],self.heymsfield_warm_mass)
            self.mass_spectrum_cold[1][td] = np.multiply(self.dsd[td],self.heymsfield_cold_mass)       
            self.mass_spectrum_conv[1][td] = np.multiply(self.dsd[td],self.heymsfield_conv_mass)       
            self.mass_spectrum_warm[0][td] = np.sum(self.mass_spectrum_warm[1][td,:])
            self.mass_spectrum_cold[0][td] = np.sum(self.mass_spectrum_cold[1][td,1:32])
            self.mass_spectrum_conv[0][td] = np.sum(self.mass_spectrum_conv[1][td,1:32])

            #Dm, Chase et al. 2020
            dm_numerator = np.sum(self.heymsfield_warm_mass*self.drop_diameter_melted*10*self.dsd[td]*self.drop_spread)
            dm_denominator = np.sum(self.heymsfield_warm_mass*self.dsd[td]*self.drop_spread)
            self.Dm[td] = np.divide(dm_numerator,dm_denominator)

            #VVD
            self.vvd[td] = np.divide(vvd_num, vvd_denom)
            #vvd_numerator = np.sum(self.heymsfield_warm_mass*self.drop_diameter_melted*self.dsd[td]*self.drop_spread)
            #vvd_denominator = np.sum(self.heymsfield_warm_mass*self.dsd[td]*self.drop_spread)
            #self.Dm[td] = np.divide(dm_numerator,dm_denominator)

            #ice mass, Matrosov 2007
            self.mass_spectrum_matrosov[1][td,0:13] = np.multiply(self.dsd[td,0:13],self.matrosov_001_02_mass)
            self.mass_spectrum_matrosov[1][td,13:30] = np.multiply(self.dsd[td,13:30],self.matrosov_02_2_mass)
            self.mass_spectrum_matrosov[1][td,30:32] = np.multiply(self.dsd[td,30:32],self.matrosov_2_mass)
            self.mass_spectrum_matrosov[0][td] = np.sum(self.mass_spectrum_matrosov[1][td,1:32])

            # IWC, Heymsfield 2005
            #mass_multiplication = np.multiply(np.divide(self.dsd[td],10^5),self.heymsfield_warm_mass)
            variable_multiplcation = np.multiply(self.dsd[td],self.heymsfield_warm_mass)
            self.IWC[1][td] = np.multiply(self.drop_spread,variable_multiplcation)
            self.IWC[1][td] = np.multiply(self.IWC[1][td],10^6)
            #self.IWC[1][td] = np.multiply(self.mass_warm[1][td],self.drop_spread)
            self.IWC[0][td] = np.sum(self.IWC[1][td,:])
            
            '''
            # LWC, Lamb and Verlinde, lecture 1 slide 41 of 510
            moments_for_LWC = np.power(np.divide(self.drop_diameter,10),3)
            coefficient = np.divide(np.pi*1,6) #density here is 1 g/cm^3 
            following_sigma = np.multiply(moments_for_LWC,self.dsd[td])
            self.LWC[1][td] = np.multiply(coefficient,following_sigma)
            self.LWC[1][td] = np.multiply(self.LWC[1][td],self.drop_spread)
            self.LWC[0][td] = np.sum(self.LWC[1][td,:])
            
            # Solving for mass using density = M/V
            moment_for_mass = np.power(np.divide(self.drop_diameter,10),3) #convert diameters from mm to cm, then put to the power of 3
            self.liquid_mass_simple[1][td,:] = (1/6)*np.pi*moment_for_mass*1 #multiply diameters^3 by pi/6 times 1 g/cm^3
            self.liquid_mass_simple[0][td] = np.sum(self.liquid_mass_simple[1][td,:])
            self.liquid_mass_spectrum_simple[1][td,:] = np.multiply(self.dsd[td],self.liquid_mass_simple[1][td]) # multiply the previous product by the respective number concentration of that diameter and time index
            self.liquid_mass_spectrum_simple[0][td] = np.sum(self.liquid_mass_spectrum_simple[1][td,:]) # sum mass spectrums from each bin to get a total number for the time index
            
            
            self.LWC[1][td] = ((1/6) * np.pi * 1e-03 * self.dsd[td] * self.drop_spread * np.array(self.drop_diameter) ** 3)
            self.LWC[0][td] = np.sum(self.LWC[1][td,:])
            '''
            # liquid mass, Williams et al. 2014
            moments_for_liquid_mass = np.power(self.drop_diameter,3) #diameters in mm to the power of 3
            coefficient = np.divide(np.pi*1,6*1000) #calculating the coefficient, density here is 1 g/cm^3 
            self.liquid_mass[1][td,:] = np.multiply(coefficient,moments_for_liquid_mass) #multiply diameters^3 by the coefficient, results in a 1x32 array, one mass for each diameter bin
            self.liquid_mass_spectrum[1][td,:] = np.multiply(self.dsd[td],self.liquid_mass[1][td]) # multiply the previous product by the respective number concentration of that diameter and time index
            self.liquid_mass_spectrum[0][td] = np.sum(self.liquid_mass_spectrum[1][td,:]) # sum mass spectrums from each bin to get a total number for the time index

            # LWC, William et al. 2014
            self.Williams_LWC[1][td,:] = np.multiply(self.liquid_mass_spectrum[1][td,:],self.drop_spread) # multiply mass spectrum by bin width to get the LWC of each bin
            self.Williams_LWC[0][td] = np.sum(self.Williams_LWC[1][td,:]) # sum LWCs from each bin to get a total LWC for the time index
            
            # equivalent reflectivity, Smith 1984
            #self.ze[0][td] = np.sum(self.ze[1][td,:])

            # #reflectivity factor (6th power of diameter, normalized for area, time)
            #self.z[1][td] = np.multiply(self.dsd[td],np.power(self.drop_diameter,6))
            #self.z[1][td] = np.multiply(self.z[1][td],self.drop_spread)
            #self.z[1][td] = np.power(self.drop_diameter,6)
            self.lwc[0][td] = np.sum(self.lwc[1][td,:])
            self.z[0][td] = np.sum(self.z[1][td,:])
            if self.z[0][td] > 0:
                self.dbz[td] = 10 * np.log10(self.z[0][td]) #z to dbz
            else:
                self.dbz[td] = float('nan')

            #self.ds.dBZ[td] = np.power(10,self.ds.dBZ[td]/10)
            #self.ds.dBZ[td] = self.ds.dBZ[td]/6
            #self.ds.dBZ[td] = 10 * np.log10(self.ds.dBZ[td])

    def saveDSD_as_nc(self,date):

        '''Takes the dsd object and saves it as a netCDF file'''
        data = xr.Dataset(data_vars=dict(
             avg_rainrate = (["time"], self.ds.rr.data), #mean Parsivel rain rate value recorded in each data chunk
             avg_dbz = (["time"], self.ds.dBZ.data), #mean Parsivel dBZ recorded in each data chunk
             #temperature = (["time"], self.proc_p2.temperature), #The mean parsivel temperature value of each data chunk
             moment_diameter = (["time"], self.dm_moments.data), #Mass-weighted mean diameter
             mass_mean_diameter = (["time"], self.Dm.data),
             
             original_dsd = (["time","diameter"], self.ds.Nd.data),
             dsd = (["time","diameter"], self.dsd.data), #The number of drops per volume of air for a given drop size interval (ie for the 32 bins)
             vvd = (["time"], self.vvd.data),

             Heymsfield_massspec_warm_bins = (["time","diameter"], self.mass_spectrum_warm[1].data),
             Heymsfield_massspec_warm_total = (["time"], self.mass_spectrum_warm[0].data),
             Heymsfield_mass_warm = (["diameter"], self.heymsfield_warm_mass.data),
             Matrosov_massspec_bins = (["time","diameter"], self.mass_spectrum_matrosov[1].data),
             Matrosov_massspec_total = (["time"], self.mass_spectrum_matrosov[0].data),
             
             IWC_bins = (["time","diameter"], self.IWC[1].data),
             IWC_total = (["time"], self.IWC[0].data),
             LWC_bins = (["time","diameter"], self.LWC[1].data),
             LWC_total = (["time"], self.LWC[0].data),
             Williams_LWC_bins = (["time","diameter"], self.Williams_LWC[1].data),
             Williams_LWC_total = (["time"], self.Williams_LWC[0].data),

             raindrop_mass_bins = (["time","diameter"], self.liquid_mass[1].data),
             raindrop_mass_spectrum_bins = (["time","diameter"], self.liquid_mass_spectrum[1].data),
             raindrop_mass_spectrum_total = (["time"], self.liquid_mass_spectrum[0].data),

             z_bins = (["time","diameter"], self.z[1].data), #Reflectivity factor from each drop size bin, value for each bin individually (plus time-dimension)
             total_z = (["time"], self.z[0].data), # reflectivity factor, total across all 32 DSD bins
             total_dbz = (["time"], self.dbz.data),
             ze_bins = (["time","diameter"], self.ze[1].data), #Reflectivity factor from each drop size bin, value for each bin individually (plus time-dimension)
             total_ze = (["time"], self.ze[0].data), # reflectivity factor, average or total across all 32 DSD bins
             #spectrum = (["time", "velocity","diameter"], np.reshape(self.ds.spectrum, (1440,32,32)))
             spectrum = (["time", "velocity","diameter"], self.ds.spectrum.data)
         ),
         coords=dict(
             time=(["time"], self.time.data), #assigning time as a coordinate
             diameter=(["diameter"], self.drop_diameter), #assigning drop diameter bins as a coordinate
             velocity=(["velocity"], self.v_parsivel), #assigning velocity bins as a coordinate
         ),
         attrs=dict(description="BUF Parsivel Data")) # adding a description for netCDF file
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
    matrosov_001_02_power = np.power(np.divide(drop_diameter[0:13],10),2)
    matrosov_02_2_power = np.power(np.divide(drop_diameter[13:30],10),2.5)
    matrosov_2_power = np.power(np.divide(drop_diameter[30:32],10),3)
    heymsfield_warm_mass = np.multiply(0.00359,heymsfield_power)
    heymsfield_cold_mass = np.multiply(0.00574,heymsfield_power)
    heymsfield_conv_mass = np.multiply(0.00630,heymsfield_power)
    matrosov_001_02_mass = np.multiply(0.003,matrosov_001_02_power)
    matrosov_02_2_mass = np.multiply(0.0067,matrosov_02_2_power)
    matrosov_2_mass = np.multiply(0.0047,matrosov_2_power)

    drop_diameter_melted = np.cbrt(np.divide(6*heymsfield_warm_mass,np.pi*1))

    
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import warnings
from datetime import date
from tqdm import tqdm
from scipy.optimize import curve_fit
import os

def calc_dsd(date):
    '''
    Funtion that wraps all three Parsivel classes together to make the
    final DSD object

    Use to get the data to make plots
    '''
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    ds = xr.open_mfdataset('/data/accp/a/mp46/Parsivel_data/'+date+'_Vt_mD_BUF.nc')

    dsd = ParsivelPSDBUF(ds)

    dsd.get_precip_params()

    dsd.saveDSD_as_nc(ds,date)

    ds.close()
    os.remove('/data/accp/a/mp46/Parsivel_data/'+date+'_Vt_mD_BUF.nc')
    return dsd

class ParsivelPSDBUF(object):

    def __init__(self,ds):
        #some data stored as tuples if applicable
        #first element = average or total across all 32 DSD bins
        #second element = value for each bin individually (plus time-dimension)
        self.ds = ds #initializing ds as the original disdrometer data file
        self.time = self.ds.time 
        self.timedim = len(ds.time) #time dimension

        self.IWC = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.Heymsfield_IWC = (np.zeros(self.timedim),np.zeros((self.timedim,32)))

        self.melted_drop_z = np.zeros(self.timedim)
        self.melted_drop_ze = np.zeros(self.timedim)

        self.Heymsfield_melted_drop_z = np.zeros(self.timedim)
        self.Heymsfield_melted_drop_ze = np.zeros(self.timedim)

        self.Manual_preciprate = (np.zeros(self.timedim),np.zeros((self.timedim,32)))
        self.Heymsfield_preciprate = (np.zeros(self.timedim),np.zeros((self.timedim,32)))

        self.Dm = np.zeros(self.timedim)
        self.Nw = np.zeros(self.timedim)
        self.mean_vt = np.zeros(self.timedim)

        self.Heymsfield_Dm = np.zeros(self.timedim)
        self.Heymsfield_Nw = np.zeros(self.timedim)

        self.Heymsfield_effective_density = np.zeros(self.timedim)
        self.effective_density = np.zeros(self.timedim)

        self.Heymsfield_SLR = np.zeros(self.timedim)
        self.SLR = np.zeros(self.timedim)

        print('Variables Initialized')

    def get_precip_params(self):
        print('Calculations in process')
        #Takes a processed_parsivel instance and calculates the precip params

        a_coefficient = np.array(self.ds.m_D_coefficient)
        a_coefficient = a_coefficient[:, np.newaxis]
        a_coefficient = np.tile(a_coefficient, (1, 32))

        # Calculating particle mass using parsivel measured diameters and the 'a' coefficient dervied by the gauge data
        # using a power law where b = 2
        diameter_squared = np.power(np.divide(np.array(self.drop_diameter),10),2)
        particle_mass = np.multiply(a_coefficient,diameter_squared)

        # Dmelted, Chase et al. 2021
        drop_diameter_melted = np.cbrt(np.divide(6*particle_mass,np.pi*1)) # 1 g/cm^3
        drop_diameter_melted_mm = drop_diameter_melted*10
        
        # Calculating the bin sizes after converting to Dmelted
        drop_diameter_meltedbinwidth1 = (drop_diameter_melted_mm[:,1]-drop_diameter_melted_mm[:,0])/2
        drop_diameter_meltedbinwidth2 = (drop_diameter_melted_mm[:,12]-drop_diameter_melted_mm[:,11])/2
        drop_diameter_meltedbinwidth3 = (drop_diameter_melted_mm[:,17]-drop_diameter_melted_mm[:,16])/2
        drop_diameter_meltedbinwidth4 = (drop_diameter_melted_mm[:,22]-drop_diameter_melted_mm[:,21])/2
        drop_diameter_meltedbinwidth5 = (drop_diameter_melted_mm[:,27]-drop_diameter_melted_mm[:,26])/2
        drop_diameter_meltedbinwidth6 = (drop_diameter_melted_mm[:,31]-drop_diameter_melted_mm[:,30])/2
        drop_diameter_meltedbinwidths = np.hstack([drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth1, drop_diameter_meltedbinwidth2, drop_diameter_meltedbinwidth2, drop_diameter_meltedbinwidth2, drop_diameter_meltedbinwidth2, drop_diameter_meltedbinwidth2, drop_diameter_meltedbinwidth3, drop_diameter_meltedbinwidth3, drop_diameter_meltedbinwidth3, drop_diameter_meltedbinwidth3, drop_diameter_meltedbinwidth3, drop_diameter_meltedbinwidth4, drop_diameter_meltedbinwidth4, drop_diameter_meltedbinwidth4, drop_diameter_meltedbinwidth4, drop_diameter_meltedbinwidth4, drop_diameter_meltedbinwidth5, drop_diameter_meltedbinwidth5, drop_diameter_meltedbinwidth5, drop_diameter_meltedbinwidth5, drop_diameter_meltedbinwidth5, drop_diameter_meltedbinwidth6, drop_diameter_meltedbinwidth6])*2
        drop_diameter_meltedbinwidths = np.reshape(drop_diameter_meltedbinwidths, (1440, 32))

        # Calculating the mass of melted snow using Dmelted and the 'a' coefficient dervied by the gauge data
        diameter_squared_melted = np.power(np.divide(drop_diameter_melted_mm,10),2.1)
        mass_melted_snow = np.multiply(a_coefficient,diameter_squared_melted)

        ################################################################################################################################################################################################
        # Same as the process above, but using the Heymsfield 2010, warm topped cloud m-D relationship
        Heymsfield_diameter_squared = np.power(np.divide(np.array(self.drop_diameter),10),2.1)
        Heymsfield_particle_mass = np.multiply(0.00359,Heymsfield_diameter_squared)

        Heymsfield_drop_diameter_melted = np.cbrt(np.divide(6*Heymsfield_particle_mass,np.pi*1)) # 1 g/cm^3
        Heymsfield_drop_diameter_melted_mm = Heymsfield_drop_diameter_melted*10
        
        Heymsfield_drop_diameter_meltedbinwidth1 = (Heymsfield_drop_diameter_melted_mm[1]-Heymsfield_drop_diameter_melted_mm[0])/2
        Heymsfield_drop_diameter_meltedbinwidth2 = (Heymsfield_drop_diameter_melted_mm[12]-Heymsfield_drop_diameter_melted_mm[11])/2
        Heymsfield_drop_diameter_meltedbinwidth3 = (Heymsfield_drop_diameter_melted_mm[17]-Heymsfield_drop_diameter_melted_mm[16])/2
        Heymsfield_drop_diameter_meltedbinwidth4 = (Heymsfield_drop_diameter_melted_mm[22]-Heymsfield_drop_diameter_melted_mm[21])/2
        Heymsfield_drop_diameter_meltedbinwidth5 = (Heymsfield_drop_diameter_melted_mm[27]-Heymsfield_drop_diameter_melted_mm[26])/2
        Heymsfield_drop_diameter_meltedbinwidth6 = (Heymsfield_drop_diameter_melted_mm[31]-Heymsfield_drop_diameter_melted_mm[30])/2
        Heymsfield_drop_diameter_meltedbinwidths = np.hstack([Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth1, Heymsfield_drop_diameter_meltedbinwidth2, Heymsfield_drop_diameter_meltedbinwidth2, Heymsfield_drop_diameter_meltedbinwidth2, Heymsfield_drop_diameter_meltedbinwidth2, Heymsfield_drop_diameter_meltedbinwidth2, Heymsfield_drop_diameter_meltedbinwidth3, Heymsfield_drop_diameter_meltedbinwidth3, Heymsfield_drop_diameter_meltedbinwidth3, Heymsfield_drop_diameter_meltedbinwidth3, Heymsfield_drop_diameter_meltedbinwidth3, Heymsfield_drop_diameter_meltedbinwidth4, Heymsfield_drop_diameter_meltedbinwidth4, Heymsfield_drop_diameter_meltedbinwidth4, Heymsfield_drop_diameter_meltedbinwidth4, Heymsfield_drop_diameter_meltedbinwidth4, Heymsfield_drop_diameter_meltedbinwidth5, Heymsfield_drop_diameter_meltedbinwidth5, Heymsfield_drop_diameter_meltedbinwidth5, Heymsfield_drop_diameter_meltedbinwidth5, Heymsfield_drop_diameter_meltedbinwidth5, Heymsfield_drop_diameter_meltedbinwidth6, Heymsfield_drop_diameter_meltedbinwidth6])*2

        Heymsfield_diameter_squared_melted = np.power(np.divide(Heymsfield_drop_diameter_melted_mm,10),2.1)
        Heymsfield_mass_melted_snow = np.multiply(0.00359,Heymsfield_diameter_squared_melted)
        
        for td in list(range(0,len(self.ds.time))): # for each time index
            print(td)
    
            # Mean of parsivel velocity measurements for each 32x32 matrix
            spectrum_sum = np.sum(np.sum(self.ds.spectrum[td, :, :],axis=1))
            self.mean_vt[td] = np.sum(np.sum(self.ds.spectrum[td, :, :],axis=1) * self.v_parsivel)/spectrum_sum
            
            # reflectivity factor (6th power of diameter, normalized for area, time)
            # you talked this one over with Steve and the values seem reasonable
            self.melted_drop_z[td] = np.sum(self.ds.dsd[td,:] * ((drop_diameter_melted_mm[td,:])**6) * drop_diameter_meltedbinwidths[td,:])
            self.melted_drop_ze[td] = (0.176 / 0.93) * self.melted_drop_z[td]

            self.Heymsfield_melted_drop_z[td] = np.sum(self.ds.dsd[td] * ((Heymsfield_drop_diameter_melted_mm)**6) * Heymsfield_drop_diameter_meltedbinwidths)
            self.Heymsfield_melted_drop_ze[td] = (0.176 / 0.93) * self.Heymsfield_melted_drop_z[td]

            #Dm, Chase et al. 2020, Eq 14
            # as per Chase et al 2020, conversion to melted diameter bin widths is not needed
            dm_numerator = np.sum(mass_melted_snow[td] * drop_diameter_melted_mm[td] * self.ds.dsd[td] * self.drop_spread)
            dm_denominator = np.sum(mass_melted_snow[td] * self.ds.dsd[td] * self.drop_spread)
            self.Dm[td] = dm_numerator / dm_denominator # units: mm

            dm_numerator = np.sum(Heymsfield_mass_melted_snow * Heymsfield_drop_diameter_melted_mm * self.ds.dsd[td] * self.drop_spread)
            dm_denominator = np.sum(Heymsfield_mass_melted_snow * self.ds.dsd[td] * self.drop_spread)
            self.Heymsfield_Dm[td] = dm_numerator / dm_denominator # units: mm

            # Liquid Equivalent Normalized Intercept Parameter, Nw, Chase et al, 2020, Eq 7
            # as per Chase et al 2020, conversion to melted diameter bin widths is not needed
            coefficient = (np.power(4,4))/6
            Nw_numerator = np.power(np.sum((np.power(drop_diameter_melted_mm[td],3) * self.ds.dsd[td] * self.drop_spread)),5)
            Nw_denominator = np.power(np.sum((np.power(drop_diameter_melted_mm[td],4) * self.ds.dsd[td] * self.drop_spread)),4)
            self.Nw[td] = coefficient*(Nw_numerator/Nw_denominator)

            Nw_numerator = np.power(np.sum((np.power(Heymsfield_drop_diameter_melted_mm,3) * self.ds.dsd[td] * self.drop_spread)),5)
            Nw_denominator = np.power(np.sum((np.power(Heymsfield_drop_diameter_melted_mm,4) * self.ds.dsd[td] * self.drop_spread)),4)
            self.Heymsfield_Nw[td] = coefficient*(Nw_numerator/Nw_denominator)
     
            # IWC, Chase et al 2021, Eq 8
            #self.IWC_Chase[0][td] = (self.Nw[td] * np.power(self.Dm[td],4) * 1 * np.pi) / 256 
            #self.IWC_Chase[0][td] = np.divide(self.IWC_Chase[0][td], 1000 ) # outputs IWC in g/m^3

            # IWC, Heymsfield 2005, Eq 2, note: the equation has a multiplier of 1e6 but through unit analysis, I deduced that's only necessary if N(D) is inputted in cgs units (so it would have units of cm^-1 cm^-3 instead of mm^-1 m^-3)
            self.IWC[1][td] = self.ds.dsd[td] * particle_mass[td] * self.drop_spread
            self.IWC[0][td] = np.sum(self.IWC[1][td, :]) # outputs IWC in g/m^3

            self.Heymsfield_IWC[1][td] = self.ds.dsd[td] * Heymsfield_particle_mass * self.drop_spread
            self.Heymsfield_IWC[0][td] = np.sum(self.Heymsfield_IWC[1][td, :]) # outputs IWC in g/m^3
            
            # Precipitation Rate, need to find a reference for this, it's in my research ntoes but I'd like a more conrete reference to literature
            self.Manual_preciprate[1][td] = self.IWC[1][td] * self.ds.terminal_velocities[td,:] * 3.6
            self.Manual_preciprate[0][td] = np.sum(self.Manual_preciprate[1][td, :]) #outputs in mm/hr

            self.Heymsfield_preciprate[1][td] = self.Heymsfield_IWC[1][td] * self.ds.terminal_velocities[td,:] * 3.6
            self.Heymsfield_preciprate[0][td] = np.sum(self.Heymsfield_preciprate[1][td, :]) #outputs in mm/hr

            # effective density, Heymsfield et al. 2004, eq 5
            V_for_timestep = self.volume_of_sphere * self.ds.dsd[td] * self.drop_spread
            self.Heymsfield_effective_density[td] = self.Heymsfield_IWC[0][td]/np.sum(V_for_timestep) # outputs in g/cm^3
            self.effective_density[td] = self.IWC[0][td]/np.sum(V_for_timestep) # outputs in g/cm^3

            # snow-to-liquid ratio (SLR), equivalent to ratio of the density of liquid water to the average density of the snow - Milbrandt et al. 2012
            self.Heymsfield_SLR[td] = self.Heymsfield_effective_density[td] / 1 # units: g/cm^3 over g/cm^3 so essentially unitless
            self.SLR[td] = self.effective_density[td] / 1 # units: g/cm^3 over g/cm^3 so essentially unitless


        print('Calculations done')

    def saveDSD_as_nc(self,ds,date):
        
        print('Saving netCDF file for '+date+'...')
        ds["mass_mean_diameter"] = (["time"], self.Dm.data)
        ds["mean_vt"] = (["time"], self.mean_vt.data)
        ds["Nw"] = (["time"], self.Nw.data)

        ds["IWC_total"] = (["time"], self.IWC[0].data)
        ds["Heymsfield_IWC_total"] = (["time"], self.Heymsfield_IWC[0].data)

        ds["Manual_preciprate_total"] = (["time"], self.Manual_preciprate[0].data)
        ds["melted_drop_ze"] = (["time"], self.melted_drop_ze.data)

        ds["Heymsfield_mass_mean_diameter"] = (["time"], self.Heymsfield_Dm.data)
        ds["Heymsfield_Nw"] = (["time"], self.Heymsfield_Nw.data)

        ds["Heymsfield_preciprate_total"] = (["time"], self.Heymsfield_preciprate[0].data)
        ds["Heymsfield_melted_drop_ze"] = (["time"], self.Heymsfield_melted_drop_ze.data)

        ds["Heymsfield_effective_density"] = (["time"], self.Heymsfield_effective_density.data)
        ds["effective_density"] = (["time"], self.effective_density.data)

        ds["Heymsfield_SLR"] = (["time"], self.Heymsfield_SLR.data)
        ds["SLR"] = (["time"], self.SLR.data)
        
        ds.to_netcdf('/data/accp/a/mp46/Parsivel_data/' + date + '_Calculated_Params_BUF.nc')  # saving netCDF to specified path
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
    
    velocities = np.reshape(np.array(v_parsivel), (32,1))

    volume_of_sphere = np.divide(np.pi,6) * np.power(np.array(drop_diameter)/10,3)
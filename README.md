# BUF Parsivel Processing Scripts

This repository includes function that handle and process Parsivel data from the UIUC SCAMP Deployment in Buffalo, NY.

In short, it filters and calculates the DSDs, fits terminal velocity relationships and derives mass-dimension relationships (from the gauge) to the data, and calculates the following parameters using the DSDs and fitted relationships (quick note, each parameter listed is calculated twice, once with the m-D relationships derived from the gauge and again using the warm-topped m-D relationship from Heymsfield et al. 2010):
1. Ice Water Content
2. Reflectivity factor using melted equivalent diameter
3. Equivalent reflectivity factor
6. Mass-weighted mean diameter 
7. Liquid Equivalent Normalized Intercept Parameter (Nw)
8. Precipitation Rate
9. Effective density
10. Snow-to-liquid ratios

Note: all literature references for each calculation can be found in the Python scripts

## Preprocessing using PyDSD

As data files get written, they automatically get processed using the PyDSD repository, particularly using [PyDSD/pydsd/io/ParsivelReader.py](https://github.com/josephhardinee/PyDSD/blob/master/pydsd/io/ParsivelReader.py). The overall function of this script is to...

1. read raw 10s parsivel files,
2. convert the time to an Epoch time
3. adds parsivel measured/outputted variables to a Python dictionary
4. Applies a data quality matrix from Ali Tokay to DSD

Each outputted dictionary then gets converted to an xarray dataset externally (not using PyDSD) and saved to a local directory as a netCDF file.

The rest of the processing flow is an adaptation from the repository created by Joe Boomgard-Zagrodnik, a postdoc a WSU, for the OLYMPEX field campaign. Joe provided the foundation I needed to create my own parsivel processing flow, which was necessary considering I knew next to nothing about python when I started grad school. 

## ParsivelPSDBUF_ND

The purpose of this script is to calculate N(D) using the methods presented in Tokay et al. 2014 and change the integration time. 

1. opens and reads the all the netCDF files for a specific date as a xarray dataset
2. resamples data from an integration time of 10s to 1min (to be consistent with gauge/wx sensor obs)
3. Initializes variables
4. calculates N(D)
5. creates a new dictionary
6. Converts dictionary to a dataset and saves as a netCDF file (1 per date)

Parameters included in the saved netCDF include:
1. N(D) (aka dsd in the script)
   - coordinates: time, diameter
2. Parsivel spectrum (32x32 matrix)
   - coordinates: time, velocity, diameter

## ParsivelPSDBUF_FittingPowerLaws_30mins

This perhaps is the more exciting part of the process because we're fitting power law relationships to the PSDs to get terminal velocity relationships and mass-dimension relationships unique to our dataset. Very briefly, one of the benefits of SCAMP is that we can use multiple instruments to retrieve microphysical information about the snow particles. One way is using the gauge data to derive m-D relationships for our dataset (m-D relationships are typically in the form of a power law, so we can use gauge data to back out a value for 'a' at each time step). This function accomplishes this task through the following workflow:

1. Ingests the previously outputted netCDF file (as dataset) for a date in addition to the gauge and wx sensor dataset.
2. Initializes variables
3. Locates the times in the gauge and wx sensor dataset that correspond to the date of the netCDF PSD.
4. For each time step, a power law relationship between Vt and D is found using scipy curvefit.
  - this power law is then applied to calculate Vt for each D at each time steps
5. A rolling mean over a 30min window is applied to the gauge data subsetted by indices found in step 2.
6. Back out 'a' using using the gauge precipitation rate, N(D), dD, Vt found in step 3, and D^2 (assuming b = 2)
7. Add the newly found power law relationships to the xarray dataset and save as a netCDF file

Parameters included in the saved netCDF include:
1. N(D) (aka dsd in the script)
   - coordinates: time, diameter
2. Parsivel spectrum (32x32 matrix)
   - coordinates: time, velocity, diameter
3. calculated terminal velocities
   - coordinates: time, diameter
4. 'a' coefficient in m-D relationship
   - coordinates: time (one a per PSD)
  
## ParsivelPSDBUF_CCalculations

Remember that long list of parameters from the introductory section of this README? Well this is where all that magic happens! All those parameters are calculated in this script and while I'll briefly go over the calculations here, references and any notes about the calculations are detailed in the script by the code.

Note: tuples are used in this script to allow the user to see contributions from each drop bin range to parameter total values (I swear this will make more sense later) where the contribution from each drop bin range is included as the second element. Parameters that are stored as tuples will be noted below.

1. Initializes variables

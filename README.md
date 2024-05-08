# BUF Parsivel Processing Scripts

This repository includes functions that handle and process Parsivel data from the UIUC SCAMP Deployment in Buffalo, NY.

In short, it filters and calculates the DSDs, fits terminal velocity relationships and derives mass-dimension relationships (from the gauge) for the data, and calculates the following parameters using the DSDs and fitted relationships (quick note, each parameter listed is calculated twice, once with the m-D relationships derived from the gauge and again using the warm-topped m-D relationship from Heymsfield et al. 2010):
1. Mean particle velocity
2. Reflectivity factor using melted equivalent diameter
3. Equivalent reflectivity factor
4. Mass-weighted mean diameter
5. Liquid equivalent normalized intercept parameter (N<sub>w</sub>)
6. Ice water content
7. Precipitation rate
8. Effective density
9. Snow-to-liquid ratios

Note: all literature references for each calculation can be found in the Python scripts and I run the scripts using the Jupyter Notebook provided.

## Preprocessing using PyDSD

As data files get written, they automatically get minorly processed using the PyDSD repository, particularly using [PyDSD/pydsd/io/ParsivelReader.py](https://github.com/josephhardinee/PyDSD/blob/master/pydsd/io/ParsivelReader.py). The overall function of this script is to:

1. Read raw 10s parsivel files
2. Convert the time to an Epoch time
3. Add parsivel measured/outputted variables to a Python dictionary

Each outputted dictionary then gets converted to an xarray dataset externally (not using PyDSD) and saved to a local directory as an hourly netCDF file.

The rest of the processing flow is an adaptation from the repository created by Joe Boomgard-Zagrodnik, a postdoc at WSU, for the OLYMPEX field campaign ([PyOLYMPEX](https://github.com/joejoezz/PyOLYMPEX.git)). 

## ParsivelPSDBUF_ND

The purpose of this script is to calculate N(D) (aka the DSD/PSD) using the methods presented in Tokay et al. 2014, and change the integration time. 

1. Opens and read the all the netCDF files for a specific date as a xarray dataset
2. Resamples data from an integration time of 10s to 1min (to be consistent with gauge/wx sensor obs)
3. Initializes variables
4. Calculates the PSD
6. Creates a new dictionary
7. Converts dictionary to a dataset and saves as a netCDF file (1 per date)

Parameters included in the saved netCDF include:
1. PSD (aka dsd in the script)
   - coordinates: time, diameter
2. Parsivel spectrum (32x32 matrix)
   - coordinates: time, velocity, diameter

## ParsivelPSD_FittingPowerLaws_30mins

This perhaps is the more exciting part of the process because we're fitting power law relationships to the PSDs to get terminal velocity relationships and mass-dimension relationships unique to our dataset. Very briefly, one of the benefits of SCAMP is that we can use multiple instruments to retrieve microphysical information about the snow particles. One way is using the gauge data to derive m-D relationships for our dataset (m-D relationships are typically in the form of a power law, so we can use gauge data to back out a value for _a_ at each time step). This function accomplishes this task through the following workflow:

1. Ingests the previously outputted netCDF file (as dataset) for a date in addition to the gauge and wx sensor dataset
2. Initializes variables
3. Locates the times in the gauge and wx sensor dataset that correspond to the date of the netCDF PSD
4. For each time step, a power law relationship between V<sub>t</sub> and D is found using scipy curvefit
   - this power law is then applied to calculate Vt for each D at each time steps
5. A rolling mean over a 30min window is applied to the gauge data subsetted by indices found in step 2
6. Back out _a_  using the gauge precipitation rate, N(D), dD, V<sub>w</sub> found in step 3, and D<sup>2</sup> (assuming b = 2)
7. Add the newly found power law relationships to the xarray dataset and save as a netCDF file

Parameters included in the saved netCDF include:
1. PSD (aka dsd in the script)
   - coordinates: time, diameter
2. Parsivel spectrum (32x32 matrix)
   - coordinates: time, velocity, diameter
3. Calculated terminal velocities
   - coordinates: time, diameter
4. _a_ coefficient in m-D relationship
   - coordinates: time (one a per PSD)
  
As previously mentioned, this function takes both gauge and Parsivel data as inputs which sometimes poses an issue. To be more descriptive, there are instances when the gauge didn't record data at a time that the Parsivel may have, and vice versa. There is only one date I have encountered this for (20221116) out of the 17 I've processed. In this case, the gauge missed 4 minutes of observations, which didn't allow my to do some array manipulation. What I did in this instance was remove the missing times from the Parsivel dataset and then proceeded to process the data. When I wrote the code to do this, I found it really slowed down the processing so in the Jupyter notebook, I switched out `ParsivelPSD_FittingPowerLaws_30mins` for `ParsivelPSD_FittingPowerLaws_AdaptiveTime` which has the code to handle issues such as this.
  
## ParsivelPSD_Calculations

Remember that long list of parameters from the introductory section of this README? Well this is where all that magic happens! All those parameters are calculated in this script and while I'll briefly go over the calculations here, references and any notes about the calculations are detailed in the script by the code.

Note: tuples are used in this script to allow the user to see contributions from each drop bin range to parameter total values where the contribution from each drop bin range is included as the second element. Parameters that are stored as tuples will be noted below and only the first element is outputted in the netCDF.

Workflow:
1. Initializes variables
2. Calculates particle mass, melted equivalent diameter, and the mass of the particles using the melted equivalent diameter (done for both the gauge derived m-D relationship and the warm-topped m-D relationship from Heymsfield et al. 2010)
3. Calculates PSD parameters (will be detailed below)
4. Add the newly calculated parameters to the xarray dataset and save as a netCDF file
   - all parameters listed above in the intro section are calculated in this script, but not all are saved/outputted (the only thing that isn't added to the dataset is Z)
6. Removes the netCDF file from the second script

Parameters included in the netCDF include:
1. Mean particle velocity
   - coordinates: time
3. Equivalent reflectivity factor
   - coordinates: time
4. Mass-weighted mean diameter
   - coordinates: time
5. Liquid equivalent normalized intercept parameter (N<sub>w</sub>)
   - coordinates: time
6. Ice water content
   - coordinates: time
   - initialized as tuple
7. Precipitation rate
   - coordinates: time
   - initialized as tuple
8. Effective density
   - coordinates: time
9. Snow-to-liquid ratios
   - coordinates: time

Note: each parameter listed (EXCEPT mean velocity) above is calculated twice, once with the m-D relationships derived from the gauge and again using the warm-topped m-D relationship from Heymsfield et al. 2010). 

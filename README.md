# BUF Parsivel Processing Scripts

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

## Preprocessing using PyDSD

As data files get written, they automatically get processed using the PyDSD repository, particularly using [PyDSD/pydsd/io/ParsivelReader.py](https://github.com/josephhardinee/PyDSD/blob/master/pydsd/io/ParsivelReader.py). The overall function of this script is to...

1. read raw 10s parsivel files,
2. convert the time to an Epoch time
3. adds parsivel measured/outputted variables to a Python dictionary
4. Applies a data quality matrix from Ali Tokay

This outputted dictionary then gets converted to an xarray dataset externally (not using PyDSD) and saved to a local directory as one netCDF file (we go from 5760 10s files to one 24 hr file, which makes it way easier to handle going forward).


## ParsivelPSDBUF_ND

After using PyDSD, 
The general flow is that Reads in daily Parsivel netCDFs and calculates n(D) using 

## ParsivelPSDBUF_FittingPowerLaws_30mins

This perhaps is the more exciting part of the process because we're fitting power law relationships to the PSDs to get terminal velocity relationships and mass-dimension relaitonships unique to our dataset.

1. Ingests the ND file for a date in addition to the gauge and wx sensor dataset.
2. Locates the times in the gauge and wx sensor dataset that correspond to the date of the PSDs.
3. For each time step,  

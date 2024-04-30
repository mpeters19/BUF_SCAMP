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

As data files get written, they automatically get processed using the PyDSD repository. While that repo offers a plethora of 
parameter calculations, it is primarily used to aggregate the raw hourly netCDFs into one daily netCDF and for it's matrix filtering 
capabilities.

Once I do a deeper dive, I will update this part.

## ParsivelPSDBUF_ND

After using PyDSD, 
The general flow is that Reads in daily Parsivel netCDFs and calculates n(D) using 

## ParsivelPSDBUF_FittingPowerLaws_30mins

This perhaps is the more exciting part of the process because we're fitting power law relationships to the PSDs to get terminal velocity relationships and mass-dimension relaitonships unique to our dataset.

1. Ingests the ND file for a date in addition to the gauge and wx sensor dataset.
2. Locates the times in the gauge and wx sensor dataset that correspond to the date of the PSDs.
3. For each time step,  

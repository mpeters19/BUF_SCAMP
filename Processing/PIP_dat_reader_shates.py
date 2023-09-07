"""
Original pip functions written to work with pip files from haukeliseter and kiruna
(opening and working with .dat files)


julia shates
10 feb 2021

"""
import numpy.ma
import numpy as np
import datetime
import netCDF4 as nc
import glob

def open_dist_data(fdatname,year,month,day,loc):
    """ Open up .dat files for PIP DSD or VVD files and save a number of variables as numpy arrays. The function also fills in missing data with -9999. (but masks it). Arrays with time dimensions have all 1440 minutes of the day. 
*~*~ usage ~*~*:
inputs: year, month, day, loc
outputs: average distribution for that day, dist per minute, center bin, edge bin, time variable, instrument time
instrument time is recorded to track what minutes were in the original file. (and VVD files have some issues with repeated time stamps)
ex:
[DIST_avg, distperminute, bin_cen, bin_edge, timevar, instrument_time] = open_dist_data(fdatname, 2018, 1, 13,'kir') 

loc options: 'kir','hauk','mqt'
this function can be used for distribution variables for velocity and particle size distributions"""
    year = int(year)
    month = int(month) 
    day = int(day) 
    DIST_avg = []
    distperminute = []
    bin_cen = []
    bin_edge = []
    timevariable = []
    instrument_time = []
    loc = loc
    start = [datetime.datetime(int(year), int(month), int(day), 0, 0)]

    dataread = []
    linelen = []
    # based on the number of columns in the line, append to a different var.
    f = open(fdatname,'r')
    for i,line in enumerate(f.readlines()):
        linelen.append(len(line.split('\t')))

    f = open(fdatname,'r')
    for i,line in enumerate(f.readlines()):
        if len(line.split('\t')) >5:
            dataread.append(line.split('\t'))
        else:
            print(line)
            #continue

    Bin_edge = dataread[1]
    Bin_edge = Bin_edge[5:]
    # make sure each term is actually a float/ not a string
    bin_edge = []
    for binpoint in Bin_edge:
        bin_edge.append(float(binpoint))
    if bin_edge[0] == 0.:
        print('removing a bin edge point')
        bin_edge = bin_edge[1:] #remove the first point only for PSD! not VVD! it's zero for psd file.

    dBin = dataread[2]
    dBin = dBin[5:]

    Bin_cen = dataread[3]
    Bin_cen= Bin_cen[5:]
    bin_cen = []
    # make sure each term is actually a float
    for binpoint in Bin_cen:
        bin_cen.append(float(binpoint))
    if float(Bin_edge[0]) == 0.:
        bin_cen = bin_cen[1:] #remove the first point...if we removed it before. 

    DIST_avg = dataread[0]
    print(DIST_avg[0:4])
    if float(Bin_edge[0]) == 0.:
        DIST_avg = DIST_avg[6:]
    else:
        DIST_avg = DIST_avg[5:]
    distavg = []
    for binpoint in DIST_avg:
        distavg.append(float(binpoint))
    DIST_avg = np.array(distavg)

    #long list thing. should be the length of the minutes for which we have data
    dist_minute = dataread[4:]

    # the time miutes that we have data
    timestamped = []
    for element in dist_minute:
        timestamped.append(float(element[0]))

    missing = np.argwhere(np.array(timestamped) == -99.)
    missing = missing[::-1] #reverse list... so indices of them don't change as we remove them

    if len(missing)>= 2:
        for missingpoint in missing:
            #remove potentially missing elements from the timevar
            del(timestamped[missingpoint])
            #remove the same potentially missing elements from the DSD_minute var as well
            del(dist_minute[missingpoint])
    else:
        #remove potentially missing elements from the timevar
        del(timestamped[int(missing)])
        #remove the same potentially missing elements from the DSD_minute var as well
        del(dist_minute[int(missing)])


# THIS PIECE WORKs. 
    #get the DSD info and ignore preamble stuff
    DIST_minute = []
    for element in dist_minute:
        if float(Bin_edge[0]) == 0.:
            interm = element[6:]
        else: 
            interm = element[5:]
        floating = []
        for elemental in interm:
            floating.append(float(elemental))
        DIST_minute.append(floating)

#BELOW HERE is less straightforward. if there are bugs in code, check this spot
# we have missing data during the day... use masking to deal 
    timealt = list(np.copy(start))
    for timestep in range(1440-1):
        timealt.append(timealt[timestep] + datetime.timedelta(minutes = 1))
    timevariable = timealt
    if not timestamped:
        print('there is no data on this day at all... but we still have a file. so here we are in this section')
        # let the first point in the day be the beginning spot
        distpermin = np.zeros((np.shape(timealt)[0], np.shape(bin_edge)[0]))
        distpermin.fill(-9999.)
        DIST_avg = np.zeros((np.shape(bin_edge)))
        DIST_avg.fill(-9999.)
    else:
        # this section deals with days that have no missing data or just a couple of missing points
        timestamped = np.array(timestamped)
        timevar = nc.num2date(timestamped[:],'days since ' + str(year)+'-' + str(month)+'-'+str(day)+' 00:00:00 UTC') 
        uptime = []
        for i,moment in enumerate(timevar):
            uptime.append(moment.replace(microsecond=0))
            if uptime[i].second == 59:
                uptime[i] = uptime[i] + datetime.timedelta(seconds = 1)

        instrument_time = uptime
        distpermin = []

        # fill in for missing data 
        zterm = np.zeros((1,np.shape(DIST_minute)[1]))
        zterm.fill(-9999.)
        for i, intertime in enumerate(timealt):
            checkpoint = np.argwhere(np.array(uptime) == intertime)
            if len(checkpoint) >1:
                print('uhoh: we have duplicate time stamps. using time variable from PSD file.')
                distpermin.extend(zterm)
                """newcheckpoint1 = handle_duplicate_time_stamps(loc, intertime)
                distpermin.append(DIST_minute[newcheckpoint1])"""
            else:
                if checkpoint:
                    distpermin.append(DIST_minute[int(checkpoint)])
                if not checkpoint:
                    distpermin.extend(zterm)
    
    distperminute = np.array(distpermin)
    distperminute = numpy.ma.masked_less_equal(distpermin,0.)
    DIST_avg = numpy.ma.masked_less_equal(DIST_avg,0.)
    return DIST_avg, distperminute, bin_cen,bin_edge, timevariable, instrument_time 

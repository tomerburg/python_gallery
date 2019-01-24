"""
NOHRSC Snow Plotting Script
Written by Tomer Burg
Last revised 1/19/2019
Python 3.6

This script retrieves 24-hour gridded snow accumulation files from NOHRSC, which
uses a blend of official COOP/COCORAHS/ASOS snow observations, radar QPE, and
model QPF, among other variables, and plots the cumulative sum of snow
accumulations over a range of dates.

For more information on NOHRSC gridded snow analyses:
https://www.nohrsc.noaa.gov/snowfall/
"""

#Import the necessary libraries
import os, sys
import urllib.request
import numpy as np
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as col

#==============================================================================
# Specify settings below
#==============================================================================

#Specify start & end dates of time range in YYYYMMDD format
#(all dates are valid at 1200 UTC)
start_date = "20170314"
end_date = "20170316"

#Specify directory to store files in. If blank, the script will use the current directory.
work_dir = ""

#Type of plot to do: (default is sum)
#"sum" = sum of snow accumulations over the time frame
#"std" = standard deviation of snow over the time frame
#"max" = maximum 24-hour snow accumulation over the time frame
plot_type = "sum"

#==============================================================================
# Step 0. Download data into directory
#==============================================================================

#Ensure working directory has a "/" at the end
if len(work_dir) > 0 and work_dir[-1] != "/": work_dir[-1] = "/"

#Convert start & end dates to datetime objects
sdate = dt.datetime.strptime(start_date+"12","%Y%m%d%H")
edate = dt.datetime.strptime(end_date+"12","%Y%m%d%H")

#Because the files cover snow accumulation for the *preceding* 24 hours, the
#start date will be incremented forward by 24 hours
sdate += dt.timedelta(hours=24)

#Empty list to store dates in, which will be used to open these files later
files = []

#Iterate through range of dates provided
while sdate <= edate:
    
    #Format current date as string
    strdate = dt.datetime.strftime(sdate,'%Y%m%d%H')
    yyyymm = dt.datetime.strftime(sdate,'%Y%m')
    
    #Format URL to download from NOHRSC
    url = f"http://www.nohrsc.noaa.gov/snowfall_v2/data/{yyyymm}/sfav2_CONUS_24h_{strdate}.nc"
    savepath = f"{work_dir}{strdate}.nc"

    #Download file
    #testfile = request.URLopener()
    #testfile.retrieve(url, savepath)
    urllib.request.urlretrieve(url, savepath)
    
    #Append filepath to list of files to open
    files.append(savepath)
    
    #Increment date by 24 hours
    sdate = sdate + dt.timedelta(hours=24)
    
    print(f"Downloaded file for {strdate}")
    
#==============================================================================
# Step 1. Load files into script
#==============================================================================

#Open files using mfdataset, concatenating by time
data = xr.open_mfdataset(files,concat_dim='time')

#Only retain snow accumulation data
snow = data.Data

#==============================================================================
# Step 2. Prepare map to plot
#==============================================================================

#==============================================================================
# Step 2. Perform calculations as requested
#==============================================================================

#Error check the provided plot_type
plot_types = ['sum','std','max']
if plot_type not in plot_types: plot_type = 'sum'

#Calculate sum
if plot_type == 'sum':
    snow_plot = snow.sum("time")

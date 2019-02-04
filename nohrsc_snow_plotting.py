"""
NOHRSC Snow Plotting Script
Written by Tomer Burg
Last revised 1/19/2019
Python 3.6

This script downloads 24-hour gridded snow accumulation files from NOHRSC, which
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
import cartopy
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from cartopy import util as cu

#==============================================================================
# Specify settings below
#==============================================================================

#Specify start & end dates of time range in YYYYMMDD format
#(all dates are valid at 1200 UTC)
start_date = "20170312"
end_date = "20170315"

start_date = "20180301"
end_date = "20180303"

#Specify directory to store files in. If blank, the script will use the current directory.
work_dir = ""

#Type of plot to do: (default is sum)
#"sum" = sum of snow accumulations over the time frame
#"std" = standard deviation of snow over the time frame
#"max" = maximum 24-hour snow accumulation over the time frame
#"days" = days of 24-hour snow accumulation above threshold
plot_type = "sum"

#If plot type is "days", specify the snow accumulation threshold desired
snow_thres = 1.0

#==============================================================================
# Step 0. Download data into directory
#==============================================================================

#Ensure working directory has a "/" at the end
if len(work_dir) > 0 and work_dir[-1] != "/": work_dir[-1] = "/"

#Convert start & end dates to datetime objects
sdate = dt.datetime.strptime(start_date+"12","%Y%m%d%H")
edate = dt.datetime.strptime(end_date+"12","%Y%m%d%H")
idate = dt.datetime.strptime(start_date+"12","%Y%m%d%H")

#Because the files cover snow accumulation for the *preceding* 24 hours, the
#start date will be incremented forward by 24 hours
idate += dt.timedelta(hours=24)

#Empty list to store dates in, which will be used to open these files later
files = []

#Iterate through range of dates provided
while idate <= edate:
    
    #Format current date as string
    strdate = dt.datetime.strftime(idate,'%Y%m%d%H')
    yyyymm = dt.datetime.strftime(idate,'%Y%m')
    
    #Format URL to download from NOHRSC
    url = f"http://www.nohrsc.noaa.gov/snowfall_v2/data/{yyyymm}/sfav2_CONUS_24h_{strdate}.nc"
    savepath = f"{work_dir}{strdate}.nc"

    #Download file, if it doesn't exist locally
    if os.path.isfile(savepath) == False:
        urllib.request.urlretrieve(url, savepath)
        print(f"---> Downloaded file for {strdate}")
        
    #Append filepath to list of files to open
    files.append(savepath)
    
    #Increment date by 24 hours
    idate +=+ dt.timedelta(hours=24)
    
#==============================================================================
# Step 1. Load files into script
#==============================================================================

#Open files using mfdataset, concatenating by time
data = xr.open_mfdataset(files,concat_dim='time')

#Only retain snow accumulation data & convert from mm to inches
snow = data.Data * 39

#Retrieve lat & lon arrays
lat = snow.lat.values
lon = snow.lon.values

print("---> Loaded files into script")

#==============================================================================
# Step 2. Prepare map to plot
#==============================================================================

#Create a Lambert conformal projection object, centered over the US
lon1 = -99.0
lat1 = 35.0
slat = 35.0

bound_n = 50.0
bound_s = 21.5
bound_w = -122.0
bound_e = -72.5
proj_lcc = ccrs.LambertConformal(central_longitude=lon1,
                             central_latitude=lat1,
                             standard_parallels=[slat])

#Create figure and axis
fig = plt.figure(figsize=(18,12),dpi=125)
ax = plt.axes(projection=proj_lcc)
ax.set_extent([bound_w,bound_e,bound_s,bound_n])

#Draw coastlines
ax.coastlines(resolution='50m', color='black', linewidths=1.2)

#Add country borders
ax.add_feature(cfeat.BORDERS, linewidths=1.2, edgecolor='black')

#Add land/lake/ocean masking
land_mask = cfeat.NaturalEarthFeature('physical', 'land', '50m',
                                    edgecolor='face',
                                    facecolor=cfeat.COLORS['land'])
sea_mask = cfeat.NaturalEarthFeature('physical', 'ocean', '50m',
                                    edgecolor='face',
                                    facecolor=cfeat.COLORS['water'])
lake_mask = cfeat.NaturalEarthFeature('physical', 'lakes', '50m',
                                    edgecolor='face',
                                    facecolor=cfeat.COLORS['water'])
ax.add_feature(sea_mask,zorder=0)
ax.add_feature(land_mask,zorder=0)
ax.add_feature(lake_mask,zorder=0)

#Add state borders
state_borders = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lakes',
                                    scale='50m', facecolor='none')
ax.add_feature(state_borders, edgecolor='black', zorder=4, linewidth=0.5)
print("---> Set up cartopy projection & geography")

#==============================================================================
# Step 2. Perform calculations and make plot
#==============================================================================

#Error check the provided plot_type
plot_types = ['sum','std','max','days']
if plot_type not in plot_types: plot_type = 'sum'

#Calculate sum
if plot_type == 'sum':
    snow_plot = snow.sum("time").values
    clevs = [0.5,1,2,3,4,6,8,12,16,20,24,30,36,48,60,72]
    cmap = col.ListedColormap(['#bdd7e7','#6baed6','#3182bd','#08519c','#082694','#ffff96','#ffc400','#ff8700','#db1400','#9e0000','#690000','#ccccff','#9f8cd8','#7c52a5','#561c72','#40dfff'])
    norm = col.BoundaryNorm(clevs,cmap.N)
    plot_title = "Cumulative Snow Accumulation (in)\n"
    
    #Locate and label maximum value in the grid
    max_val = np.nanmax(snow_plot)
    idx = np.where(snow_plot==max_val)
    max_lon = lon[idx[1]]
    max_lat = lat[idx[0]]
    plt.plot(max_lon,max_lat,'*',ms=22,mec='k',mew=2.7,fillstyle='none',transform=ccrs.PlateCarree(),label='Maximum Snow:\n%0.1f"'%(max_val))

#Maximum 24-hour snowfall
if plot_type == 'max':
    snow_plot = snow.max("time").values
    clevs = [0.5,1,2,3,4,6,8,12,16,20,24,30,36,48,60,72]
    cmap = col.ListedColormap(['#bdd7e7','#6baed6','#3182bd','#08519c','#082694','#ffff96','#ffc400','#ff8700','#db1400','#9e0000','#690000','#ccccff','#9f8cd8','#7c52a5','#561c72','#40dfff'])
    norm = col.BoundaryNorm(clevs,cmap.N)
    plot_title = "Maximum 24-Hour Snow Accumulation (in)\n"
    
    #Locate and label maximum value in the grid
    max_val = np.nanmax(snow_plot)
    idx = np.where(snow_plot==max_val)
    max_lon = lon[idx[1]]
    max_lat = lat[idx[0]]
    plt.plot(max_lon,max_lat,'*',ms=22,mec='k',mew=2.7,fillstyle='none',transform=ccrs.PlateCarree(),label='Maximum Snow:\n%0.1f"'%(max_val))
    
#Standard deviation of snowfall, for all cases where snow over 0.1" fell
if plot_type == 'std':
    
    #replace snow below 0.1 with nan's
    snow = snow.values
    snow[snow < 0.1] = np.nan
    snow_plot = np.nanstd(snow,axis=0)
    
    clevs = np.arange(1,np.nanmax(snow_plot),1)
    cmap = plt.cm.YlGnBu
    norm = col.BoundaryNorm(clevs,cmap.N)
    plot_title = "Standard Deviation of 24-Hour Snow Accumulation (in)\n"
    
#Days of snow above specified threshold
if plot_type == 'days':
    
    #replace snow below threshold with nan's
    snow_plot = snow.values
    snow_plot[snow_plot < snow_thres] = 0.0
    snow_plot[snow_plot >= snow_thres] = 1.0
    snow_plot = np.nansum(snow_plot,axis=0)
    
    clevs = np.arange(1,np.nanmax(snow_plot)+1,1)
    cmap = plt.cm.YlGnBu
    norm = col.BoundaryNorm(clevs,cmap.N)
    plot_title = f"Days of 24-Hour Snowfall Above {snow_thres} inches\n"
    
    #Contour 1 day threshold
    ax.contour(lon,lat,snow_plot,[1],colors='k',linewidths=0.2,transform=ccrs.PlateCarree())
    
    #Locate and label maximum value in the grid
    max_val = np.nanmax(snow_plot)
    idx = np.where(snow_plot==max_val)
    max_lon = lon[idx[1]]
    max_lat = lat[idx[0]]
    plt.plot(max_lon,max_lat,'*',ms=22,mec='k',mew=2.7,fillstyle='none',transform=ccrs.PlateCarree(),label='Maximum Days:%s'%(str(int(max_val))))


print("---> Finished calculating data for plot")

#Add legend
plt.legend(loc=3, prop={'size': 14})

#Add date range to plot title
plot_title += dt.datetime.strftime(sdate,'%Y-%m-%d %H00 UTC') + " - " + dt.datetime.strftime(edate,'%Y-%m-%d %H00 UTC')

#Fill contour the snowfall data
cs = plt.contourf(lon,lat,snow_plot,clevs,cmap=cmap,norm=norm,transform=ccrs.PlateCarree())
plt.colorbar(cs,shrink=0.79,pad=0.01,ticks=clevs)

#Plot title
plt.title(plot_title,fontsize=16,fontweight='bold')

#Save image and close
fname = dt.datetime.strftime(sdate,'%Y%m%d') + "_" + dt.datetime.strftime(edate,'%Y%m%d') + "_" + plot_type + ".png"
plt.savefig(fname,bbox_inches='tight')
plt.close()
data.close()
print("Done!")

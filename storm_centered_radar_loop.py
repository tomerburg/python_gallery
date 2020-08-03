"""
Storm-Centered MRMS Radar Loop Script
Written by Tomer Burg
Last revised 8/2/2020
Python 3.7

This script accesses MRMS mosaic base reflectivity data from UCAR's thredds
server, provided it is within the last month, and generates storm-centered
images every 10 minutes.

UCAR Thredds server (used to retrieve MRMS data):
https://thredds.ucar.edu/thredds/catalog.html

Tropycal python package (used to retrieve tropical cyclone data):
https://tropycal.github.io/tropycal/

Cartopy wrapper (used to generate maps):
https://github.com/tomerburg/Map
"""

#Import the necessary libraries
import os, sys
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
import cartopy.crs as ccrs
import tropycal.tracks as tracks
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.dates as mdates
import matplotlib.patheffects as path_effects

#Import cartopy wrapper
from Map import Map

#==============================================================================
# Specify settings below
#==============================================================================

#Specify requested storm from Tropycal (storm name, storm year)
storm_tuple = ('isaias',2020)

#Loop start time (if commented out, default is start of the storm)
start_date = dt.datetime(2020,8,1,19)

#Loop end time (if commented out, default is end of the storm)
end_date = dt.datetime(2020,8,1,20)

#Output directory to save images in (default is current directory)
output_directory = ""

#If True, images will be exactly storm-centered. If False, images will be slightly offset
#to the north of the storm center.
storm_center = True

#==============================================================================
# Define general functions to be used in this script
#==============================================================================

def temporal_interpolation(values, orig_times, target_times, kind='quadratic'):
    r"""
    Interpolate an array of values between two time arrays.
    
    Parameters
    ----------
    values : 1D array
        1D array of values to be interpolated (e.g., wind, lat, lon)
    orig_times : 1D array
        1D array of datetime objects corresponding to values
    target_times : 1D array
        1D array of datetime objects corresponding to the desired output
    kind : str, optional
        Type of interpolation. Default is quadratic. Can be linear.
    
    Returns
    -------
    1D array of values corresponding to the target time array.
    """
    
    #Interpolate array and get output values
    f = interp.interp1d(mdates.date2num(orig_times),values,kind=kind)
    ynew = f(mdates.date2num(target_times))
    
    #Return new array
    return ynew

def find_nearest(array,value):
    r"""
    Find the index of the closest value in an array to the desired value.
    
    Parameters
    ----------
    array : 1D array
        1D array of values.
    value : int or float
        Integer or float to search for in the array.
    
    Returns
    -------
    int
        Integer corresponding to the index of the closest value in array to value.
    """
    
    return np.abs(array - value).argmin()

def storm_cat(storm_wind):
    r"""
    Determine the Saffir-Simpson Hurricane Wind Scale (SSHWS) category based on sustained wind.
    
    Parameters
    ----------
    storm_wind : int or float
        Sustained wind in knots.
    
    Returns
    -------
    str
        String corresponding to the storm category.
    """
    
    if storm_wind < 34: return "Tropical Depression"
    if storm_wind < 64: return "Tropical Storm"
    if storm_wind < 83: return "Category 1 Hurricane"
    if storm_wind < 96: return "Category 2 Hurricane"
    if storm_wind < 113: return "Category 3 Hurricane"
    if storm_wind < 136: return "Category 4 Hurricane"
    return "Category 5 Hurricane"

#==============================================================================
# Define colormap related functions to be used in this script
#==============================================================================

def add_alpha(r1,g1,b1,a):
    r"""
    Add opacity to an RGB value.
    """
    r2 = 255
    g2 = 255
    b2 = 255

    r3 = r2 + (r1-r2)*a
    g3 = g2 + (g1-g2)*a
    b3 = b2 + (b1-b2)*a
    
    return [r3,g3,b3]

def rgb(r,g,b):
    r"""
    Return a hex string corresponding to RGB values.
    """
    r = int(r)
    g = int(g)
    b = int(b)
    return '#%02x%02x%02x' % (r, g, b)
    
def getColor(val,rng,col1,col2):
    r"""
    Get a color within a linear range.
    """
    r1,g1,b1 = col1
    r2,g2,b2 = col2
    
    #168,0,168 to 0,0,230 = -168,0,+62
    #rng = (-16 - value)
    rdif = float(r2 - r1)
    gdif = float(g2 - g1)
    bdif = float(b2 - b1)
    
    r3 = r2 + (-1.0 * val * (rdif / float(rng)))
    g3 = g2 + (-1.0 * val * (gdif / float(rng)))
    b3 = b2 + (-1.0 * val * (bdif / float(rng)))

    return rgb(r3,g3,b3)

def reflectivity(clevs):
    r"""
    Return a list of colors for reflectivity (dBZ).
    """
    
    r1,g1,b1 = add_alpha(160,168,180,0.05)
    r2,g2,b2 = add_alpha(67,94,159,0.85)
    
    colors = []
    for value in clevs:
        
        if value < 0:
            colors.append(rgb(r1,g1,b1))
            #colors.append(rgb(120,137,174))
            
        elif value < 12:
            rng = 12
            val = (12 - value)
            colors.append(getColor(val,rng,[r1,g1,b1],[r2,g2,b2]))
            
        elif value < 20:
            rng = 8
            val = (20 - value)
            colors.append(getColor(val,rng,[r2,g2,b2],[111,214,232]))
            
        elif value < 24:
            rng = 4
            val = (24 - value)
            colors.append(getColor(val,rng,[111,214,232],[17,213,24]))
            
        elif value < 33: #37
            rng = 9
            val = (33 - value)
            colors.append(getColor(val,rng,[17,213,24],[9,94,9]))
            
        elif value < 40: #42.5
            rng = 7
            val = (40 - value)
            colors.append(getColor(val,rng,[9,94,9],[255,226,0]))
            
        elif value < 50:
            rng = 10
            val = (50 - value)
            colors.append(getColor(val,rng,[255,226,0],[255,128,0]))
            
        elif value < 60:
            rng = 10
            val = (60 - value)
            colors.append(getColor(val,rng,[255,0,0],[113,0,0]))
            
        elif value < 65:
            rng = 5
            val = (65 - value)
            colors.append(getColor(val,rng,[255,245,255],[255,146,255]))
            
        elif value < 70:
            rng = 5
            val = (70 - value)
            colors.append(getColor(val,rng,[255,117,255],[225,11,227]))
            
        elif value < 75:
            rng = 5
            val = (75 - value)
            colors.append(getColor(val,rng,[178,0,255],[99,0,214]))

        else:
            colors.append(rgb(99,0,214))

    return colors


#==============================================================================
# Retrieve storm data from Tropycal
#==============================================================================

#Retrieve HURDATv2 dataset, including Best Track for current year data
basin = tracks.TrackDataset(include_btk=True)

#Retrieve data for the requested storm
storm = basin.get_storm(storm_tuple)

#==============================================================================
# Interpolate HURDATv2 or best track data to every 10 minutes
#==============================================================================

#Determine loop start and end dates
try:
    start_date
except:
    start_date = storm.date[0]
    
try:
    end_date
except:
    end_date = storm.date[-1]

#Construct array of target times using MRMS data
times_target = []
while start_date <= end_date:
    times_target.append(start_date)
    start_date += dt.timedelta(minutes=10)

#Construct data arrays interpolated to every 10 minutes
lats = temporal_interpolation(storm.lat,storm.date,times_target)
lons = temporal_interpolation(storm.lon,storm.date,times_target)
vmax = temporal_interpolation(storm.vmax,storm.date,times_target,kind='linear')
mslp = temporal_interpolation(storm.mslp,storm.date,times_target,kind='linear')

#==============================================================================
# Interpolate HURDATv2 or best track data to every 10 minutes
#==============================================================================

#Iterate over all time steps
for idx,(lon,lat,i_vmax,i_mslp,i_time) in enumerate(zip(lons,lats,vmax,mslp,times_target)):

    #Try reading in data. If there's an error, it will continue the script
    try:

        #Read reflectivity data from UCAR's thredds server using xarray
        url = f"https://thredds.ucar.edu/thredds/dodsC/grib/NCEP/MRMS/BaseRef/MRMS_BaseReflectivity_{i_time.strftime('%Y%m%d_%H')}00.grib2"
        ds = xr.open_dataset(url)

        #Identify closest time in this grib2 file
        times = ds.time.values
        times_num = np.array([mdates.date2num(i) for i in ds.time.values])
        closest_idx = find_nearest(times_num,mdates.date2num(i_time))
        new_time = times[closest_idx]
        new_time = pd.to_datetime(new_time)

        #Retrieve reflectivity data, subset to storm-centered coordinates
        refl = ds.MergedBaseReflectivityQC_altitude_above_msl.sel(lat=slice(lat+4.5,lat-4.5),lon=slice(lon+360.0-6.5,lon+360.0+6.5))
        refl = refl.isel(altitude_above_msl=0,time=closest_idx)
        refl_lat = refl.lat.values
        refl_lon = refl.lon.values

        #Close dataset
        ds.close()

        #----------------------------------------------------------------------

        #Create storm-centered map projection
        m = Map('LambertConformal',
                central_longitude = lon,
                central_latitude = lat,
                standard_parallels = [lat],
                res = 'h',
        )

        #Create figure
        fig = plt.figure(figsize=(14,9))
        ax = plt.axes(projection=m.proj)

        #Draw continents and oceans
        ax.background_patch.set_facecolor('#EDFBFF')
        m.filllakes('#EDFBFF',zorder=2)
        m.fillcontinents('#FBF5EA',zorder=1)

        #Draw geography
        m.drawcoastlines(zorder=4)
        m.drawcountries(zorder=4)
        m.drawstates(zorder=4)

        #----------------------------------------------------------------------

        #Contour fill reflectivity data
        clevs = np.arange(5,75,1.0)
        clevsbar = np.arange(5,75,5.0)
        cmap = col.ListedColormap(reflectivity(clevs))
        norm = col.BoundaryNorm(clevs,cmap.N)

        #Replace all values less than 5 dBZ with NaN
        refl = refl.values
        refl = np.nan_to_num(refl)
        refl[refl < 5.0] = np.nan

        #Contour fill plot and add colorbar
        cs = ax.pcolormesh(refl_lon-360.0,refl_lat,refl,vmin=clevs[0],vmax=clevs[-1],cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),zorder=3)
        m.colorbar(cs,ticks=clevsbar)

        #----------------------------------------------------------------------

        #Set map extent to storm-centered
        if storm_center == True:
            ax.set_extent((-500*1000.0, 500*1000.0, -325*1000.0, 325*1000.0), crs=m.proj)
        else:
            ax.set_extent((-500*1000.0, 500*1000.0, -250*1000.0, 400*1000.0), crs=m.proj)

        #Plot title
        plt.title("MRMS Base Reflectivity (dBZ)",loc='left',fontsize=16,fontweight='bold')
        plt.title(f"{new_time.strftime('%H%M UTC %d %B %Y')}",loc='right',fontsize=14)

        #Add current intensity
        new_text = f"{storm_cat(i_vmax)}\nMax Wind: {str(int(i_vmax))} kt\nMin MSLP: {str(int(i_mslp))} hPa"
        plt.text(0.02,0.97,new_text,ha='left',va='top',transform=ax.transAxes,fontsize=14,
                 bbox=dict(facecolor='#ffffff', edgecolor='black', boxstyle='round,pad=0.3', alpha=0.5))
        new_text = storm_cat(i_vmax)
        a = plt.text(0.02,0.97,new_text,ha='left',va='top',transform=ax.transAxes,fontsize=14)
        a.set_path_effects([path_effects.Stroke(linewidth=0.4, foreground='k'),path_effects.Normal()])

        #Show plot and close
        file_path = os.path.join(output_directory,f"{idx}.png")
        plt.savefig(file_path,bbox_inches='tight')
        plt.close()

    #If there's an error, simply continue to the next time step
    except:

        print(f"Error for time {new_time}")
        continue

print("Done!")

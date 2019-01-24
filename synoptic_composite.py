"""
CFSR Synoptic Composite Plotting Script
Written by Tomer Burg
Last revised 1/24/2019
Python 3.6

This script retrieves CFSR data from UAlbany's THREDDS server, and plots
a composite of various synoptically-relevant variables using MetPy and Cartopy with
a standard mercator projection. This script plots data for a single time step, although
it is possible to edit this into a loop over multiple time steps.
"""

#Import the necessary libraries
import numpy as np
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as col
import scipy.ndimage as ndimage
from scipy.ndimage.filters import minimum_filter, maximum_filter
import matplotlib.patheffects as path_effects
import matplotlib.gridspec as gridspec

import cartopy
from cartopy import crs as ccrs
import cartopy.feature as cfeature
from cartopy import util as cu

import metpy
import metpy.calc as calc

#========================================================================================================
# Specify settings below
#========================================================================================================

#Enter the date and hour to visualize data for, in YYYYMMDDHH format (e.g., 3/13/1993 18 UTC = 1993031318)
plot_date = "1993031318"

#Choose the lat/lon range over which to get data for and plot the map.
map_subset = 1 #if zero, then the following 4 variables will be ignored and a global plot will be made.
lat_south = 20.0
lat_north = 60.0
lon_west = 230.0-360.0
lon_east = 300.0-360.0

#========================================================================================================
# Step 1. Load the variables from UAlbany's opendap server
#========================================================================================================

#To make parsing the date easier, convert it into a datetime object and get it into various formats
date_obj = dt.datetime.strptime(plot_date,'%Y%m%d%H')
yyyy = date_obj.year
mm = date_obj.month
dd = date_obj.day
yyyymmdd = dt.datetime.strftime(date_obj,'%Y%m%d')

#For loading multiple files in
files = []

#Loop through each variable specified
variables = ['u','v','g','t','pmsl','pwat']
for var in variables:
    
    #Append file into list of files to open
    filepath = "http://thredds.atmos.albany.edu:8080/thredds/dodsC/CFSR/%s/%s.%s.0p5.anl.nc"%(yyyy,var,yyyy)
    files.append(filepath)
    
#Load in the variable(s) as an xarray Dataset and assign them into "data"
data = xr.open_mfdataset(files)

print("Loaded dataset from UAlbany thredds server")

#========================================================================================================
# Step 2. Subset the data by pressure level, time and geographic region
#========================================================================================================

#Subset by the time
data = data.sel(time=date_obj)

#Subset data by latitude and longitude, if requested
if map_subset == 1:
    data = data.sel(lat=slice(lat_south-10,lat_north+10),lon=slice(lon_west-10,lon_east+10))
    
#Get lat and lon arrays for this dataset:
lat = data.lat.values
lon = data.lon.values

print("Subsetted data")

#========================================================================================================
# Step 3. Create a Basemap plotting figure and add geography
#========================================================================================================

#Create a Plate Carree projection object
proj_ccrs = ccrs.Miller(central_longitude=0.0)

#Create figure and axes for main plot and colorbars
fig = plt.figure(figsize=(18,12),dpi=125)
gs = gridspec.GridSpec(12, 36, figure=fig) #[ytop:ybot, xleft:xright]
ax = plt.subplot(gs[:, :-1],projection=proj_ccrs) #main plot
ax.set_xticklabels([])
ax.set_yticklabels([])
ax2 = plt.subplot(gs[:4, -1]) #top plot
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax3 = plt.subplot(gs[4:8, -1]) #bottom plot
ax3.set_xticklabels([])
ax3.set_yticklabels([])
ax4 = plt.subplot(gs[8:, -1]) #bottom plot
ax4.set_xticklabels([])
ax4.set_yticklabels([])

#Add political boundaries and coastlines
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidths=1.2)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidths=1.2)
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidths=0.5)

#Add land/lake/ocean masking
land_mask = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                    edgecolor='face',
                                    facecolor='#e6e6e6')
sea_mask = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                    edgecolor='face',
                                    facecolor='#ffffff')
lake_mask = cfeature.NaturalEarthFeature('physical', 'lakes', '50m',
                                    edgecolor='face',
                                    facecolor='#ffffff')
ax.add_feature(sea_mask,zorder=0)
ax.add_feature(land_mask,zorder=0)
ax.add_feature(lake_mask,zorder=0)

print("Set up cartopy projection & geography")

#========================================================================================================
# Supplementary functions
#========================================================================================================

#Creates a color gradient based on passed ranges of values
class gradient:
    
    def __init__(self,*args):
        self.args = args
        
        #Set threshold levels
        self.thres = []
        self.thres_min = []
        
        #Error check arguments
        error = 0
        for arg in args:
            
            #Ensure each argument & subargument has a length of 2
            if len(arg) != 2: error = 1
            if len(arg[0]) != 2: error = 1
            if len(arg[1]) != 2: error = 1
            
            #Ensure the 2nd element of each argument is a number
            if isinstance(arg[0][1], (int, float)) == False: error = 2
            if isinstance(arg[1][1], (int, float)) == False: error = 2
            
            #Ensure that the 1st element of each argument is either a hex str or rgb tuple
            if isinstance(arg[0][0], (str, tuple)) == False: error = 3
            if isinstance(arg[1][0], (str, tuple)) == False: error = 3
            
            #Ensure gradient values are continuous
            if len(self.thres) > 0 and self.thres[-1] != arg[0][1]: error = 4
            
            #Append threshold levels
            self.thres.append(arg[0][1])
            self.thres.append(arg[1][1])
            self.thres_min.append(arg[0][1])
            
        #Ensure values are either constantly increasing or decreasing
        check_thres = np.array(self.thres)
        diff = check_thres[1:] - check_thres[:-1]
        if np.min(diff) == 0 and np.max(diff) > 0:
            pass
        elif np.min(diff) < 0 and np.max(diff) == 0:
            self.thres = self.thres[::-1]
        else:
            error = 4
        
        #Output error messages
        if error == 1: raise RuntimeError('Each argument must have 2 elements, e.g., [["#00FFFF",25.0],["#0000FF",29.0]]')
        if error == 2: raise RuntimeError('The second element must be a number, e.g., [["#00FFFF",25.0]')
        if error == 3: raise RuntimeError('The first element must be a hex string or an rgb tuple, e.g., [["#00FFFF",25.0]')
        if error == 4: raise RuntimeError('Values assigned to the gradient must be continuous, either increasing or decreasing.')
        
    #Returns the hex string corresponding to the passed rgb values
    def rgb(self,r,g,b):
        r = int(r)
        g = int(g)
        b = int(b)
        return '#%02x%02x%02x' % (r, g, b)
        
    #Computes a hex value matching up with the current position relative to the range of colors.
    #position = current position within the range of colors (e.g., 1)
    #rng = range of colors (e.g. 5, so 1/5 would be 20% of the range)
    #col1 = Starting RGB color for the range (e.g. [0,255,255])
    #col2 = Ending RGB color for the range (e.g. [0,0,255])
    def getColor(self,position,rng,col1,col2):
        
        #Retrieve r,g,b values from tuple
        r1,g1,b1 = col1
        r2,g2,b2 = col2
    
        #Get difference in each r,g,b value between start & end 
        rdif = float(r2 - r1)
        gdif = float(g2 - g1)
        bdif = float(b2 - b1)
        
        #Calculate r,g,b values for the specified position within the range
        r3 = r2 + (-1.0 * position * (rdif / float(rng)))
        g3 = g2 + (-1.0 * position * (gdif / float(rng)))
        b3 = b2 + (-1.0 * position * (bdif / float(rng)))
    
        #Return in hex string format
        return self.rgb(r3,g3,b3)

    #Finds the nearest gradient range to use
    def find_nearest(self,arr,val):
        for ival in arr[::-1]:
            if ival <= val:
                return arr.index(ival)
        
    #Create a color map based on passed levels
    def get_cmap(self,levels):
        
        #Add empty color list
        self.colors = []
        
        #Iterate through levels
        for lev in levels:
            
            #Check if level is outside of range
            if lev < self.thres[0]:
                start_hex = self.args[0][0][0]
                if "#" not in start_hex: start_hex = self.rgb(start_hex[0],start_hex[1],start_hex[2])
                self.colors.append(start_hex)
            
            elif lev > self.thres[-1]:
                end_hex = self.args[-1][1][0]
                if "#" not in end_hex: end_hex = self.rgb(end_hex[0],end_hex[1],end_hex[2])
                self.colors.append(end_hex)
                
            else:
                
                #Find closest lower threshold
                idx = self.find_nearest(self.thres_min,lev)
                
                #Retrieve start & end values
                start_value = self.args[idx][0][1]
                end_value = self.args[idx][1][1]
                
                #Calculate start and end RGB tuples, if passed as hex
                start_hex = self.args[idx][1][0]
                end_hex = self.args[idx][0][0]
                if "#" in start_hex:
                    start_hex = start_hex.lstrip('#')
                    end_hex = end_hex.lstrip('#')
                    start_rgb = tuple(int(start_hex[i:i+2], 16) for i in (0, 2 ,4))
                    end_rgb = tuple(int(end_hex[i:i+2], 16) for i in (0, 2 ,4))
                else:
                    start_rgb = start_hex
                    end_rgb = end_hex
    
                #Get hex value for the color at this point in the range
                nrange_color = (end_value - start_value)
                idx = lev - start_value
                hex_val = self.getColor(idx,nrange_color,start_rgb,end_rgb)
                
                #Append color to list
                self.colors.append(hex_val)
        
        #Convert to a colormap and return
        self.cmap = col.ListedColormap(self.colors)
        return self.cmap
            
#Spatially smooth a 2D variable
def smooth(prod,sig):
    
    #Check if variable is an xarray dataarray
    try:
        lats = prod.lat.values
        lons = prod.lon.values
        prod = ndimage.gaussian_filter(prod,sigma=sig,order=0)
        prod = xr.DataArray(prod, coords=[lats, lons], dims=['lat', 'lon'])
    except:
        prod = ndimage.gaussian_filter(prod,sigma=sig,order=0)
    
    return prod

#Label MSLP extrema
def extrema(mat,mode='wrap',window=50):
        mn = minimum_filter(mat, size=window, mode=mode)
        mx = maximum_filter(mat, size=window, mode=mode)
        return np.nonzero(mat == mn), np.nonzero(mat == mx)

def mslp_label(mslp,lat,lon):
    
    #Determine an appropriate window given the lat/lon grid resolution
    res = lat[1] - lat[0]
    nwindow = int(9.5 / res)
    mslp = np.ma.masked_invalid(mslp)
    local_min, local_max = extrema(mslp, mode='wrap', window=nwindow)
    
    #Determine axis boundaries
    xmin, xmax, ymin, ymax = ax.get_extent()
    lons2d, lats2d = np.meshgrid(lon, lat)
    transformed = proj_ccrs.transform_points(proj_ccrs, lons2d, lats2d)
    x = transformed[..., 0]
    y = transformed[..., 1]
    
    #Get location of extrema on grid
    xlows = x[local_min]; xhighs = x[local_max]
    ylows = y[local_min]; yhighs = y[local_max]
    lowvals = mslp[local_min]; highvals = mslp[local_max]
    yoffset = 0.022*(ymax-ymin)
    dmin = yoffset
    
    #Plot low pressures
    xyplotted = []
    for x,y,p in zip(xlows, ylows, lowvals):
        if x < xmax-yoffset and x > xmin+yoffset and y < ymax-yoffset and y > ymin+yoffset:
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin: #,fontweight='bold'
                a = ax.text(x,y,'L',fontsize=28,
                        ha='center',va='center',color='r',fontweight='normal')
                b = ax.text(x,y-yoffset,repr(int(p)),fontsize=14,
                        ha='center',va='top',color='r',fontweight='normal')
                a.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='black'),
                       path_effects.SimpleLineShadow(),path_effects.Normal()])
                b.set_path_effects([path_effects.Stroke(linewidth=1.0, foreground='black'),
                       path_effects.SimpleLineShadow(),path_effects.Normal()])
                xyplotted.append((x,y))
                
    #Plot high pressures
    xyplotted = []
    for x,y,p in zip(xhighs, yhighs, highvals):
        if x < xmax-yoffset and x > xmin+yoffset and y < ymax-yoffset and y > ymin+yoffset:
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                a = ax.text(x,y,'H',fontsize=28,
                        ha='center',va='center',color='b',fontweight='normal')
                b = ax.text(x,y-yoffset,repr(int(p)),fontsize=14,
                        ha='center',va='top',color='b',fontweight='normal')
                a.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='black'),
                       path_effects.SimpleLineShadow(),path_effects.Normal()])
                b.set_path_effects([path_effects.Stroke(linewidth=1.0, foreground='black'),
                       path_effects.SimpleLineShadow(),path_effects.Normal()])
                xyplotted.append((x,y))

#========================================================================================================
# Step 4. Fill contours
#========================================================================================================

#--------------------------------------------------------------------------------------------------------
# 850-hPa temperature
#--------------------------------------------------------------------------------------------------------

#Get the data for this variable
temp = data['t'].sel(lev=850)
temp.metpy.convert_units('degC')

#Specify contour settings
clevs = np.arange(-40,40,1)
cmap = plt.cm.jet
extend = "both"

#Contour fill this variable
norm = col.BoundaryNorm(clevs,cmap.N)
cs = ax.contourf(lon,lat,temp,clevs,cmap=cmap,norm=norm,extend=extend,transform=proj_ccrs,alpha=0.1)

print("Filled contours for 850-hPa temperature")

#--------------------------------------------------------------------------------------------------------
# PWAT
#--------------------------------------------------------------------------------------------------------

#Get the data for this variable
pwat = data['pwat']
pwat.metpy.convert_units('mm')

#Specify contour settings
clevs = np.arange(20,71,0.5)

#Define a color gradient for PWAT
pwat_colors = gradient([[(255,255,255),0.0],[(255,255,255),20.0]],
               [[(205,255,205),20.0],[(0,255,0),34.0]],
               [[(0,255,0),34.0],[(0,115,0),67.0]])
cmap = pwat_colors.get_cmap(clevs)
extend = "max"

#Contour fill this variable
norm = col.BoundaryNorm(clevs,cmap.N)
cs = ax.contourf(lon,lat,pwat,clevs,cmap=cmap,norm=norm,extend=extend,transform=proj_ccrs,alpha=0.9)

#Add a color bar
cbar = plt.colorbar(cs,cax=ax2,shrink=0.75,pad=0.01,ticks=[20,30,40,50,60,70])

print("Filled contours for PWAT")

#--------------------------------------------------------------------------------------------------------
# 250-hPa wind
#--------------------------------------------------------------------------------------------------------

#Get the data for this variable
u = data['u'].sel(lev=250)
v = data['v'].sel(lev=250)
wind = calc.wind_speed(u,v)

#Specify contour settings
clevs = [40,50,60,70,80,90,100,110]
cmap = col.ListedColormap(['#99E3FB','#47B6FB','#0F77F7','#AC97F5','#A267F4','#9126F5','#E118F3','#E118F3'])
extend = "max"

#Contour fill this variable
norm = col.BoundaryNorm(clevs,cmap.N)
cs = ax.contourf(lon,lat,wind,clevs,cmap=cmap,norm=norm,extend=extend,transform=proj_ccrs)

#Add a color bar
cbar = plt.colorbar(cs,cax=ax3,shrink=0.75,pad=0.01,ticks=clevs)

print("Filled contours for 250-hPa wind")

#--------------------------------------------------------------------------------------------------------
# 500-hPa smoothed vorticity
#--------------------------------------------------------------------------------------------------------

#Get the data for this variable
u = data['u'].sel(lev=500)
v = data['v'].sel(lev=500)
dx,dy = calc.lat_lon_grid_deltas(lon,lat)
vort = calc.vorticity(u,v,dx,dy)
smooth_vort = smooth(vort, 5.0) * 10**5

#Specify contour settings
clevs = np.arange(2,20,1)
cmap = plt.cm.autumn_r
extend = "max"

#Contour fill this variable
norm = col.BoundaryNorm(clevs,cmap.N)
cs = ax.contourf(lon,lat,smooth_vort,clevs,cmap=cmap,norm=norm,extend=extend,transform=proj_ccrs,alpha=0.3)

#Add a color bar
cbar = plt.colorbar(cs,cax=ax4,shrink=0.75,pad=0.01,ticks=clevs[::2])

print("Filled contours for 500-hPa vorticity")
        
#========================================================================================================
# Step 5. Contours
#========================================================================================================

#--------------------------------------------------------------------------------------------------------
# MSLP
#--------------------------------------------------------------------------------------------------------

#Get the data for this variable
mslp = data['pmsl']
mslp = mslp / 100.0

#Specify contour settings
clevs = np.arange(960,1040+4,4)
style = 'solid' #Plot solid lines
color = 'red' #Plot lines as gray
width = 0.8 #Width of contours 0.25

#Contour this variable
cs = ax.contour(lon,lat,mslp,clevs,colors=color,linewidths=width,linestyles=style,transform=proj_ccrs,alpha=0.9)

#Include value labels
ax.clabel(cs, inline=1, fontsize=9, fmt='%d')

print("Contours complete for MSLP")

#--------------------------------------------------------------------------------------------------------
# Geopotential heights
#--------------------------------------------------------------------------------------------------------

#Get the data for this variable
hght = data['g'].sel(lev=500)
hght = hght / 10.0

#Specify contour settings
clevs = np.arange(480,612,12)
style = 'solid' #Plot solid lines
color = 'black' #Plot lines as gray
width = 2.0 #Width of contours

#Contour this variable
cs = ax.contour(lon,lat,hght,clevs,colors=color,linewidths=width,linestyles=style,transform=proj_ccrs)

#Include value labels
ax.clabel(cs, inline=1, fontsize=12, fmt='%d')

print("Contours complete for 500-hPa heights")


#--------------------------------------------------------------------------------------------------------
# Surface barbs
#--------------------------------------------------------------------------------------------------------

#Get data for wind barbs
u = data['u'].sel(lev=850)
v = data['v'].sel(lev=850)

#Plot wind barbs
quivers = ax.quiver(lon, lat, u.values, v.values, transform=proj_ccrs, regrid_shape=(38,30), scale=820, alpha=0.5)

print("Barbs complete for 850-hPa wind")

#--------------------------------------------------------------------------------------------------------
# Label highs & lows
#--------------------------------------------------------------------------------------------------------

#Label highs and lows
mslp_label(mslp,lat,lon)

#========================================================================================================
# Step 6. Add legend, plot title, then save image and close
#========================================================================================================

#Add custom legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='#00A123', lw=5),
                Line2D([0], [0], color='#0F77F7', lw=5),
                Line2D([0], [0], color='#FFC000', lw=5),
                Line2D([0], [0], color='k', lw=2),
                Line2D([0], [0], color='k', lw=0.1, marker=r'$\rightarrow$', ms=20),
                Line2D([0], [0], color='r', lw=0.8),
               ]

ax.legend(custom_lines, ['PWAT (mm)', '250-hPa Wind (m/s)', '500-hPa Vorticity', '500-hPa Height (dam)', '850-hPa Wind (m/s)', 'MSLP (hPa)'], loc=2, prop={'size':12})

#Format plot title
title = "Synoptic Composite \nValid: " + dt.datetime.strftime(date_obj,'%Y-%m-%d %H%M UTC')
st = plt.suptitle(title,fontweight='bold',fontsize=16)
st.set_y(0.92)

#Save figure and close
plt.savefig(f"synoptic_composite_{plot_date}.png",bbox_inches='tight')
plt.close()

print("Done!")

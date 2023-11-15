"""
directory: name of directory (e.g., "2m_temperature", "2pvu_theta_wind", "500mb_geopotential")
grid_spacing: "0p25" or "0p5"
"""

import os, sys
import xarray as xr
import numpy as np
import datetime as dt
import calendar

#=======================================================================================================================

directory = sys.argv[1]
grid_spacing = sys.argv[2]

#=======================================================================================================================

# Year range for climatology
years = range(1991,2021)

# Open datasets (each file is 2.9 GB, node memory constraint is 128 GB RAM)
ds = xr.open_mfdataset([f'ERA5/{grid_spacing}/{directory}/{year}.nc' for year in years])

# Dict mapping variable name to output file names for special cases
special_cases = {
    '250mb_wind': {
        'u': '250mb_u',
        'v': '250mb_v',
    },
    '500mb_wind': {
        'u': '500mb_u',
        'v': '500mb_v',
    },
    '850mb_wind': {
        'u': '850mb_u',
        'v': '850mb_v',
    },
    '10m_wind': {
        'u10': '10m_u',
        'v10': '10m_v',
    },
    '2pvu_theta_wind': {
        'pt': '2pvu_theta',
        'u': '2pvu_u',
        'v': '2pvu_v',
    },
}

# Iterate over all variables in file
variables = [var for var in ds.variables if var not in ['time', 'latitude', 'longitude']]
for variable in variables:
    print(f'Calculating climatology for: {variable}')

    # Start and end times for climatology
    start_time = dt.datetime(2020,1,1,0)
    end_time = dt.datetime(2020,12,31,18)

    # List to store all times and data arrays
    all_times = []
    all_data_mean = []
    all_data_std = []

    # Iterate over all 6-hourly time steps
    while start_time <= end_time:
        print(start_time)

        # Create a list of times to correspond to a centered 61-day window, at 24 hour time increments, and
        # weights from 0 to 1 to correspond to a triangular weight centered at the middle of the 61-day window
        times = []
        weights = []
        for year in years:

            # Current handling for February 29th is to center the time on February 28th
            if start_time.strftime('%m%d') == '0229' and not calendar.isleap(year):
                iter_time = (start_time-dt.timedelta(hours=24)).replace(year=year)
                times += [(iter_time - dt.timedelta(hours=i*24)).replace(year=year) for i in range(-30,31)]
            else:
                iter_time = (start_time).replace(year=year)
                times += [(iter_time - dt.timedelta(hours=i*24)).replace(year=year) for i in range(-30,31)]
            weights += [(-1*abs(i)+31)/31 for i in range(-30,31)]
        
        # Subset these times out of the full dataset
        ds_subset = ds.sel(time=times)

        # Calculate weighted mean and standard deviation
        weights = xr.DataArray(weights, dims='time')
        weighted_data = ds_subset[variable].weighted(weights)
        all_data_mean.append(weighted_data.mean('time').compute())
        all_data_std.append(weighted_data.std('time').compute())

        # Store time and increment
        all_times.append(start_time)
        start_time += dt.timedelta(hours=6)

    # Convert data to xarray Dataset
    da_mean = xr.DataArray(all_data_mean,coords=[all_times,ds.latitude.values,ds.longitude.values],dims=['time','latitude','longitude'])
    da_std = xr.DataArray(all_data_std,coords=[all_times,ds.latitude.values,ds.longitude.values],dims=['time','latitude','longitude'])
    ds = xr.Dataset({
        variable: da_mean,
        f'{variable}_std': da_std
    })

    # Save as netCDF file
    if directory in special_cases:
        save_name = special_cases.get(directory)[variable]
    else:
        save_name = directory + ''
    ds.to_netcdf(f'climatology/{grid_spacing}/{save_name}.nc')

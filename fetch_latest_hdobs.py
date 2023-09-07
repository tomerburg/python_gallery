from tropycal import recon

# ================================================================================
# Modify this section with your requested parameters
# ================================================================================

# Enter requested basin ("north_atlantic" or "east_pacific")
basin = 'north_atlantic'

# Enter aircraft ("usaf" or "noaa")
aircraft = 'usaf'

# ================================================================================
# Main code below
# ================================================================================

# Function to nicely format recon data
def format_recon_hdobs(data):
    print(f'Aircraft: {data["mission_id"]}')
    print(f'Observation Window: {data["start_time"].strftime("%H%M UTC")} - {data["end_time"].strftime("%H%M UTC %d %b %Y")}')
    print(f'\nMinimum MSLP (ESTIMATED): {data["min_mslp"]} hPa')
    print(f'Maximum SFMR Wind: {data["max_sfmr"]} knots')
    print(f'Maximum 30s Wind: {data["max_wspd"]} knots')
    print(f'\nMaximum Temperature: {data["max_temp"]} C')
    print(f'Maximum Dewpoint: {data["max_dwpt"]} C')

# Fetch realtime recon object
obj = recon.RealtimeRecon()

# Fetch parsed HDOBs data
hdobs_dict = obj.get_hdobs_realtime(basin=basin, aircraft=aircraft)

# Nicely format and output HDOBs data
format_recon_hdobs(hdobs_dict)

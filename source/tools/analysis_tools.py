"""
Functions to analyse the TLEs and ephemerides
"""
import os
import re
import pandas as pd
import numpy as np
import datetime
import json
from multiprocessing import Pool

#local imports
from .tletools import read_TLEs, TLE_time, sgp4_prop_TLE, combine_TLE2eph
from .conversions import jd_to_utc, utc_jd_date, midnight_jd_date, HCL_diff, dist_3d, alt_series, ecef_to_lla, eci2ecef_astropy

rmse = lambda x: np.sqrt(np.mean(np.square(x)))

def master_sgp4_ephemeris(start_date, stop_date, master_TLE, update = True, dt = 15*60, ):
    """Given a list of TLEs, a start day, and a stop day, return a list of SGP4 propagated ephemerides for each TLE in the list

    Args:
        start_day (string): start day in the format [YYYY,MM,DD, HH, MM, SS]
        stop_day (string): stop day in the format [YYYY,MM,DD, HH, MM, SS]
        TLE_list (list): list of TLEs
        dt (float): time step in seconds (default is 15 minutes)
        update (bool): if True, the ephemeris will be updated with each new TLE (default is True)

    Returns:
        ephemeris, orbit_ages: ephemeris in list format, list of orbit ages in hours that correspond to the ephemeris
    """
        #convert the start and stop dates to jd time stamps (midnight on the day)
    jd_start = utc_jd_date(start_date[2], start_date[1], start_date[0], 0, 0, 0) #midnight on the start day
    jd_stop = utc_jd_date(stop_date[2], stop_date[1], stop_date[0], 0, 0, 0) #midnight on the stop day


    # Check the date of the first and last TLEs in each list
    # If the dates are outside the range of the start and stop dates set the start and stop dates to the first and last TLE dates

    for sc in range(0, len(master_TLE), 1): # for each spacecraft's list of TLEs
        first_SC_TLEs = master_TLE[sc][0] # get the first TLE in the listfor this particular sc
        last_SC_TLEs = master_TLE[sc][-1] # get the last TLE in the list for this particular sc

        if jd_start < TLE_time(first_SC_TLEs):
            print('The specified start date is earlier than the date of the first TLE in the list')
            print('Setting the start date to the date of the first TLE in the list')
            #set the start date to the date of the first TLE in the list + 15 minutes
            jd_current_time = TLE_time(first_SC_TLEs)
            #convert to datetime object
            jd_dt = datetime.datetime(1858, 11, 17) + datetime.timedelta(jd_current_time - 2400000.5)
            # convert to day, month, year
            jd_start = midnight_jd_date(jd_dt.day, jd_dt.month, jd_dt.year)
            print('Start date set to: ', jd_start)
            print('equivalent to midnight on(dd,mm,yyy)',jd_dt.day,"-",jd_dt.month,"-",jd_dt.year)

        if jd_stop > TLE_time(last_SC_TLEs):
            print('The stop time is later than the date of the last TLE in the list')
            print('Setting the stop date to the date of the last TLE in the list')
            #set the stop date to the date of the last TLE in the list - 15 minutes
            jd_current_time = TLE_time(last_SC_TLEs)
            #convert to datetime object
            jd_dt = datetime.datetime(1858, 11, 17) + datetime.timedelta(jd_current_time - 2400000.5)
            # convert to day, month, year
            jd_stop = midnight_jd_date(jd_dt.day, jd_dt.month, jd_dt.year)
            print('Stop date set to: ', jd_start)
            print('equivalent to midnight on(dd,mm,yyy)',jd_dt.day,"-",jd_dt.month,"-",jd_dt.year)
        
        #if the start date is later than the stop date, return an error
        if jd_start > jd_stop:
            print('The start date is later than the stop date')
            print('Please check the start and stop dates')
            return
    
    if update == False: # this does not update the ephermeris with each new TLE
        #make an empty array of the same shape as master_TLE called master_ephemeris
        master_ephemeris = np.empty_like(master_TLE, dtype=object)
        
        ## One loop that just propagates all the TLEs within the given time period
        for sc in range(0, len(master_TLE), 1): # for each spacecraft's list of TLEs
            SC_TLEs = master_TLE[sc] # list of TLEs for a single spacecraft
            
            for tle in range(0, len(SC_TLEs), 1): # for each TLE in the list of TLEs for a single spacecraft
                single_TLE = SC_TLEs[tle]
            prop = sgp4_prop_TLE(single_TLE, jd_start, jd_stop, dt) # propagate the TLE over the selected time period 
            master_ephemeris[sc][tle] = prop #fill in the relevant spot in the master ephemeris array
        
    elif update == True: # this updates the ephemeris with each new TLE
        ## Another loop that combines the ephemerides for the most recent TLEs into a single ephemeris for each spacecraft
        master_ephemeris = np.zeros(len(master_TLE), dtype=object)
        orbit_ages = np.zeros(len(master_TLE), dtype=object)
    
        for sc in range(0, len(master_TLE), 1): # for each spacecraft's list of TLEs
            SC_TLEs = master_TLE[sc] # list of TLEs for a single spacecraft
            SC_eph, SC_orbit_ages = combine_TLE2eph(SC_TLEs, jd_start, jd_stop) # combine the TLEs into a single ephemeris over the determined time period
            master_ephemeris[sc] = SC_eph
            orbit_ages[sc] = SC_orbit_ages

    return master_ephemeris, orbit_ages

def master_sgp4_ephemeris_optimized(start_date, stop_date, master_TLE, update = True, dt = 15*60):
    jd_start = utc_jd_date(start_date[2], start_date[1], start_date[0], 0, 0, 0) 
    jd_stop = utc_jd_date(stop_date[2], stop_date[1], stop_date[0], 0, 0, 0) 

    # Compute TLE times once and store them
    TLE_times = [[TLE_time(TLE) for TLE in SC_TLEs] for SC_TLEs in master_TLE]

    for sc, SC_TLEs in enumerate(master_TLE):
        first_SC_TLEs_time = TLE_times[sc][0] 
        last_SC_TLEs_time = TLE_times[sc][-1] 

        if jd_start < first_SC_TLEs_time:
            jd_current_time = first_SC_TLEs_time
            jd_dt = datetime.datetime(1858, 11, 17) + datetime.timedelta(jd_current_time - 2400000.5)
            jd_start = midnight_jd_date(jd_dt.day, jd_dt.month, jd_dt.year)

        if jd_stop > last_SC_TLEs_time:
            jd_current_time = last_SC_TLEs_time
            jd_dt = datetime.datetime(1858, 11, 17) + datetime.timedelta(jd_current_time - 2400000.5)
            jd_stop = midnight_jd_date(jd_dt.day, jd_dt.month, jd_dt.year)
        
    if jd_start > jd_stop:
        print('The start date is later than the stop date')
        print('Please check the start and stop dates')
        return

    master_ephemeris = np.empty_like(master_TLE, dtype=object) if not update else np.zeros(len(master_TLE), dtype=object)
    orbit_ages = np.zeros(len(master_TLE), dtype=object)

    for sc, SC_TLEs in enumerate(master_TLE): 
        if update:
            master_ephemeris[sc], orbit_ages[sc] = combine_TLE2eph(SC_TLEs, jd_start, jd_stop)
        else:
            for tle, single_TLE in enumerate(SC_TLEs):
                master_ephemeris[sc][tle], _ = sgp4_prop_TLE(single_TLE, jd_start, jd_stop, dt)

    return master_ephemeris, orbit_ages

def TLE_pair_analyse(pair_TLE_list):
    
    """Takes an input two lists of TLEs. For the common time period over which there are TLEs in the list, the TLEs of each list will be propagated using the SGP4 propagator.
       A state vector [Julian Day, (x-eci, y-eci, z-eci),(u-eci, v-eci, w-eci)] will be output every 15 minutes for each spacecraft from midnight beggining on the first day of the available common time period of the provided TLE lists.
       These state vectors will be appended to a list for each spacecraft, and the pair of ephemerides will are returned in a list called "master_ephs".
        master_ephs,eph_alts, h_diffs, c_diffs, l_diffs, cart_pos_diffs, times

    Args:
        pair_TLE_list (list): a list of two lists containing TLEs for the two spacecraft to be compared. TLE list can be generated by the read_TLEs() function. 
        plot (bool, optional): Will plot the outputs below into various subplots. Defaults to False.
        savepath (str, optional): Path in which to save said plots. Defaults to './TLE_analysis_results/'.

    Returns:
        master_ephs: the list of ephemerides
        eph_alts: a list of altitudes for each spacecraft at each step
        h_diffs: a list of height differences between each spacecraft at each step
        c_diffs: a list of cross-track differences between each spacecraft at each step
        l_diffs: a list of along-track differences between each spacecraft at each step
        cart_pos_diffs: a list of 3D cartesian position differences between each spacecraft at each step
        times: a list of times for each step 
    """
    tlelist1 = pair_TLE_list[0]
    tlelist2 = pair_TLE_list[1]

    #find the time period which is common to both TLE lists
    start = max(TLE_time(tlelist1[0]), TLE_time(tlelist2[0]))
    end = min(TLE_time(tlelist1[-1]), TLE_time(tlelist2[-1]))

    #add one day to the start time and remove one day from the end time
    newstart = start + 1
    newend = end - 1

    start_year = (int(str(jd_to_utc(newstart))[0:4]))
    start_month = (int(str(jd_to_utc(newstart))[5:7]))
    start_day = (int(str(jd_to_utc(newstart))[8:10]))

    end_year = (int(str(jd_to_utc(newend))[0:4]))
    end_month = (int(str(jd_to_utc(newend))[5:7]))
    end_day = (int(str(jd_to_utc(newend))[8:10]))

    #generate master ephemeris for the time period
    master_ephs, orbit_ages = master_sgp4_ephemeris([start_year, start_month, start_day, 0, 0, 0], [end_year,end_month, end_day, 0, 0, 0], pair_TLE_list)

    # Extract position, times and velocities into arrays
    positions = np.zeros_like(master_ephs, dtype=object) # array of the same shape as master_ephs but filled with zeros
    pos_vels = np.zeros_like(master_ephs, dtype=object) # array of the same shape as master_ephs but filled with zeros

    #array containing the positons of the spacecraft at each time step
    for i in range (0, len(master_ephs), 1):
        pos_vecs = []
        for j in range(0, len(master_ephs[i]), 1):
            pos = master_ephs[i][j][1]
            pos_vecs.append(pos)
        positions[i] = pos_vecs
    
    #3D distance between each time step
    cart_pos_diffs = []
    for i in range(0, len(positions[0]), 1):
        diff = dist_3d(positions[0][i],positions[1][i])
        cart_pos_diffs.append(diff)

    # array of time stamps (for plotting)
    times = []
    for i in range(0, len(master_ephs[0]), 1):
        time = master_ephs[0][i][0]
        times.append(time)

    #from each of the ephemerides, extract the position vectors and the velocity vectors and append them to one array
    for i in range (0, len(master_ephs), 1):
        pos_vel_vecs = []
        for j in range(0, len(master_ephs[i]), 1):
            pos_vel = master_ephs[i][j][1:3]
            posvel_group = []
            for k in range(0, len(pos_vel), 1):
                for l in range (0, len(pos_vel[k]), 1):
                    posvel_group.append(pos_vel[k][l])
            pos_vel_vecs.append(posvel_group)
        pos_vels[i] = pos_vel_vecs

    # Project into HCL 
    h_diffs, c_diffs, l_diffs = HCL_diff(pos_vels[0], pos_vels[1])

    #altitude time series for each ephemeris in the amster ephemeris
    eph_alts = [alt_series(ephemeris) for ephemeris in master_ephs]
        
    return master_ephs,eph_alts, h_diffs, c_diffs, l_diffs, cart_pos_diffs, times, orbit_ages

def analyze_files(args):
    NORAD_file, SUP_file, NORAD_TLE_folder, SUP_TLE_folder, analysis_output_path = args
    #extract the numerical values from NORAD_file and SUP_file
    NORAD_id = re.findall(r'\d+', NORAD_file)
    SUP_id = re.findall(r'\d+', SUP_file) 
    NORAD_NORAD = NORAD_id
    SUP_NORAD = SUP_id
    if NORAD_NORAD == SUP_NORAD:
        print('A NORAD/SUP TLE pair was found for NORAD ID:', NORAD_NORAD[0])
        filename = str(NORAD_NORAD[0]) + '.csv'
        total_out_path = analysis_output_path + filename
        #set the paths to the files
        NORAD_TLE_path = NORAD_TLE_folder + NORAD_file
        SUP_TLE_path = SUP_TLE_folder + SUP_file
        #Read the TLEs
        print('Reading SUP TLEs')
        sup_read = read_TLEs(SUP_TLE_path)
        print('Reading NORAD TLEs')
        NORAD_read = read_TLEs(NORAD_TLE_path)
        #combine the TLEs into a list 
        sup_NORAD_pair = [sup_read, NORAD_read]
        #Analyse the differences 
        master_ephs, eph_alts, h_diffs, c_diffs, l_diffs, cart_pos_diffs, times, orbit_ages = TLE_pair_analyse(sup_NORAD_pair)
        # Make pandas dataframe and save into csv
        df = pd.DataFrame({'h_diffs': h_diffs, 'c_diffs': c_diffs, 'l_diffs': l_diffs, 'cart_pos_diffs': cart_pos_diffs, 'times': times, 'eph_alts_sup': eph_alts[0], 'eph_alts_norad': eph_alts[1], 'master_ephs_sup': master_ephs[0], 'master_ephs_norad': master_ephs[1], 'orbit_ages_sup': orbit_ages[0], 'orbit_ages_norad': orbit_ages[1]})
        df.to_csv(total_out_path)

def NORAD_vs_SUP_TLE_analysis(NORADS = [], analysis_output_path = 'output/TLE_analysis/'):
    NORAD_TLE_folder = 'external/NORAD_TLEs/'
    SUP_TLE_folder = 'external/SUP_TLEs/'

    NORAD_TLE_files = os.listdir(NORAD_TLE_folder)
    SUP_TLE_files = os.listdir(SUP_TLE_folder)

    if NORADS == []:
        raise ValueError('Please specify a list of NORAD IDs to analyze')
    else:
        NORAD_TLE_files = [file for file in NORAD_TLE_files if file.startswith(tuple(NORADS))]
        SUP_TLE_files = [file for file in SUP_TLE_files if file.startswith(tuple(['sup_' + NORAD for NORAD in NORADS]))]

        # Create a pool of workers
        with Pool() as pool:
            # Use starmap to apply the analyze_files function to each pair of NORAD and SUP files
            pool.map(analyze_files, [(NORAD_file, SUP_file, NORAD_TLE_folder, SUP_TLE_folder, analysis_output_path) for NORAD_file in NORAD_TLE_files for SUP_file in SUP_TLE_files])

def process_ephemeris_data(eph_str):
    """Convert the ephemeris strings from the TLE analysis files into a list of floats.

    Args:
        eph_str (str): The ephemeris string from the TLE analysis files

    Returns:
        list: The ephemeris string converted to a list of floats
    """
    eph_str = eph_str.split(' ')
    eph_time = eph_str[0][1:-1]
    pos1 = eph_str[1][1:-1]
    pos2 = eph_str[2][0:-1]
    pos3 = eph_str[3][0:-2]
    vel1 = eph_str[4][1:-1]
    vel2 = eph_str[5][0:-1]
    vel3 = eph_str[6][0:-2]
    eph_time = float(eph_time)
    pos1 = float(pos1)
    pos2 = float(pos2)
    pos3 = float(pos3)
    vel1 = float(vel1)
    vel2 = float(vel2)
    vel3 = float(vel3)
    return [eph_time, pos1, pos2, pos3, vel1, vel2, vel3]

import numpy as np
import json

def add_latlon_to_dfs(df):
    """Add latitude and longitude columns to the dataframe.
       These are calculated from the ephemeris data in the dataframe.

    Args:
        df (pandas dataframe): The dataframe to add the latitude and longitude columns to.

    Returns:
        df: The dataframe with the latitude and longitude columns added.
    """
    eci_positions = np.array(df['master_ephs_sup'].to_list())[:, 1:4]
    eci_velocities = np.array(df['master_ephs_sup'].to_list())[:, 4:7]
    mjd_times = df['times'] - 2400000.5
    ecef_positions, _ = eci2ecef_astropy(eci_positions, eci_velocities, mjd_times)
    lla_coords = np.array([ecef_to_lla(x, y, z) for x, y, z in ecef_positions])
    lats = lla_coords[:, 0]
    lons = lla_coords[:, 1]
    df['lat'] = lats
    df['lon'] = lons
    return df

def add_launch_numbers_to_df(df):
    """
    Add the launch numbers to the dataframe based on the NORAD IDs of the satellites in the dataframe.

    Args:
        df (pandas dataframe): The dataframe to add the launch numbers to.

    Returns:
        df: The dataframe with the launch numbers added.
    """
    with open('external/selected_satellites.json', 'r') as f:
        selected_satellites = json.load(f)

    NORAD_IDs = df['NORAD_ID'].to_list()
    launch_numbers = []

    for NORAD_ID in NORAD_IDs:
        launch_number = re.search(r'L(\d+)', NORAD_ID)
        if launch_number:
            launch_numbers.append(int(launch_number.group(1)))
        else:
            # Handle cases where the launch number is not found
            launch_numbers.append(None)

    df['launch'] = launch_numbers
    return df


def TLE_analysis_to_df():
    TLE_analysis_path = "output/TLE_analysis/"
    oneweb_NORAD_IDs = set(open('external/Constellation_NORAD_IDs/oneweb_NORAD_IDs.txt', 'r').read().splitlines())
    starlink_NORAD_IDs = set(open('external/Constellation_NORAD_IDs/starlink_NORAD_IDs.txt', 'r').read().splitlines())
    oneweb_dfs = []
    starlink_dfs = []

    for file in os.listdir(TLE_analysis_path):
        if file.endswith('.csv'):
            norad_id = file[:-4]

            if norad_id in oneweb_NORAD_IDs:
                df = pd.read_csv(TLE_analysis_path + file)
                df['NORAD_ID'] = norad_id
                df['constellation'] = 'OneWeb'
                oneweb_dfs.append(df)

            elif norad_id in starlink_NORAD_IDs:
                df = pd.read_csv(TLE_analysis_path + file)
                df['NORAD_ID'] = norad_id
                df['constellation'] = 'Starlink'
                starlink_dfs.append(df)

    #Process the ephemeris data if it hasn't already been processed
    for df in oneweb_dfs:
        df['master_ephs_sup'] = df['master_ephs_sup'].apply(process_ephemeris_data)
        df = add_latlon_to_dfs(df)
        df = add_launch_numbers_to_df(df)

    for df in starlink_dfs:
        df['master_ephs_sup'] = df['master_ephs_sup'].apply(process_ephemeris_data)
        df = add_latlon_to_dfs(df)
        df = add_launch_numbers_to_df(df)

    return oneweb_dfs, starlink_dfs
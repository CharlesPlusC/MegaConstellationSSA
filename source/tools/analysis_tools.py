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
from typing import List, Dict, Tuple, Optional
from astropy.time import Time
import scipy as sp
import scipy.fftpack

#local imports
from .tletools import twoLE_parse, read_TLEs, TLE_time, sgp4_prop_TLE, combine_TLE2eph, load_satellite_lists, read_spacex_ephemeris, spacex_ephem_to_dataframe, tle_convert
from .conversions import car2kep, kep2car, jd_to_utc, utc_jd_date, midnight_jd_date, HCL_diff, dist_3d, alt_series, ecef_to_lla, eci2ecef_astropy, eci2latlon, TEME_to_MEME, parse_spacex_datetime_stamps, yyyy_mm_dd_hh_mm_ss_to_jd

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

    #lat, lon only calculated for the first ephemeris (supplemental TLE data)
    # get the positions form pos_vels[0]. For each item in pos_vels[0], get the first 3 items (position) and the last next 3 items (velocity)
    positions = [pos_vels[0][i][:3] for i in range(len(pos_vels[0]))]
    velocities = [pos_vels[0][i][3:6] for i in range(len(pos_vels[0]))]

    # subtract 2400000.5 from times to convert from julian date to mjd time
    times = np.array(times)
    mjd_times = [times[i] - 2400000.5 for i in range(len(times))]
    lats, lons = eci2latlon(eci_positions=positions, eci_velocities=velocities, mjd_times=mjd_times)

    return master_ephs,eph_alts, h_diffs, c_diffs, l_diffs, cart_pos_diffs, times, orbit_ages, lats, lons

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
        master_ephs, eph_alts, h_diffs, c_diffs, l_diffs, cart_pos_diffs, times, orbit_ages, lats, lons = TLE_pair_analyse(sup_NORAD_pair)
        # Make pandas dataframe and save into csv
        df = pd.DataFrame({'h_diffs': h_diffs, 'c_diffs': c_diffs, 'l_diffs': l_diffs, 'cart_pos_diffs': cart_pos_diffs, 'times': times, 'eph_alts_sup': eph_alts[0], 'eph_alts_norad': eph_alts[1], 'master_ephs_sup': master_ephs[0], 'master_ephs_norad': master_ephs[1], 'orbit_ages_sup': orbit_ages[0], 'orbit_ages_norad': orbit_ages[1], 'lats': lats, 'lons': lons})
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

    # Build a dictionary that maps directly from NORAD_ID to launch_number
    norad_to_launch = {}
    for constellation in selected_satellites:
        for launch in selected_satellites[constellation]:
            for norad_id in selected_satellites[constellation][launch]:
                norad_to_launch[norad_id] = int(launch[1:])  # Exclude the "L" and convert to integer

    # Now, you can simply use this dictionary to add the launch numbers to the dataframe
    df['NORAD_ID'] = df['NORAD_ID'].astype(int)
    df['launch_no'] = df['NORAD_ID'].map(norad_to_launch)

    # Check if any NORAD ID was not found in the dictionary
    if df['launch_no'].isna().any():
        missing_norad_ids = df.loc[df['launch_no'].isna(), 'NORAD_ID'].tolist()
        raise ValueError(f"NORAD ID(s) {missing_norad_ids} not found in selected_satellites")

    return df
    
def process_TLE_analysis_file(file: str, TLE_analysis_path: str, oneweb_NORAD_IDs: set, starlink_NORAD_IDs: set) -> pd.DataFrame:
    # read in the TLE analysis .csv file and apply minor processing changes

    print('Reading in file: ' + file)
    norad_id = file[:-4]
    df = pd.read_csv(TLE_analysis_path + file)
    df['NORAD_ID'] = norad_id

    if norad_id in oneweb_NORAD_IDs:
        df['constellation'] = 'OneWeb'
    elif norad_id in starlink_NORAD_IDs:
        df['constellation'] = 'Starlink'

    df['master_ephs_sup'] = df['master_ephs_sup'].apply(process_ephemeris_data)
    df['UTC'] = Time(df['times'].values, format='jd', scale='utc').datetime # convert the times to UTC
    df['MJD'] = df['times'] - 2400000.5 # convert the times to MJD
    df = add_launch_numbers_to_df(df)

    # go through the dataframes in OneWeb_dfs and check if the NORAD ID is in the L6 list
    # if it is remove all the rows that have times before 2459500 because before this date they are orbit raising
    satellite_lists = load_satellite_lists()
    norad_OW_L6 = satellite_lists["oneweb"]["L6"]
    if int(df['NORAD_ID'][0]) in norad_OW_L6:
        mask = df['master_ephs_sup'].apply(lambda x: x[0] >= 2459500)
        df = df[mask]
        df.reset_index(drop=True, inplace=True)
    
    return df

def TLE_analysis_to_df(NORAD_IDs: list = None):
    """
    Analyze TLE (Two-Line Element Set) data and convert it to pandas DataFrame.
    
    This function reads TLE analysis files and categorizes them based on the constellation they belong to.
    Each file is read into a DataFrame, which is then processed and appended to the respective list.
    The function returns two lists of DataFrames, each corresponding to a specific constellation.
    
    Parameters:
    NORAD_IDs (list, optional): List of NORAD IDs to be analyzed. If not specified, all TLE analysis files are analyzed.

    Returns:
    tuple: Two lists of pandas DataFrames. The first list corresponds to the OneWeb constellation, 
    and the second to the Starlink constellation.
    """

    TLE_analysis_path = "output/TLE_analysis/"
    print("specified NORAD IDs: " + str(NORAD_IDs))

    oneweb_NORAD_IDs = set(open('external/Constellation_NORAD_IDs/oneweb_NORAD_IDs.txt', 'r').read().splitlines())
    starlink_NORAD_IDs = set(open('external/Constellation_NORAD_IDs/starlink_NORAD_IDs.txt', 'r').read().splitlines())
    oneweb_dfs, starlink_dfs = [], []

    if NORAD_IDs:
        print("Reading TLE analysis files for NORAD IDs: " + str(NORAD_IDs))
        oneweb_NORAD_IDs = [x for x in NORAD_IDs if x in oneweb_NORAD_IDs]
        starlink_NORAD_IDs = [x for x in NORAD_IDs if x in starlink_NORAD_IDs]
        not_found = set(NORAD_IDs) - (set(oneweb_NORAD_IDs) | set(starlink_NORAD_IDs))  # NORAD IDs not found in either list
        if not_found:
            raise ValueError("NORAD IDs not found: " + str(not_found))
    else:
        print("no NORAD IDs specified- analyzing all TLE analysis files")
        oneweb_NORAD_IDs = set(oneweb_NORAD_IDs)
        starlink_NORAD_IDs = set(starlink_NORAD_IDs)

        if len (oneweb_NORAD_IDs) != 0:
            print("OneWeb NORAD IDs: ", oneweb_NORAD_IDs)
        elif len (starlink_NORAD_IDs) != 0:
            print("Starlink NORAD IDs: " ,starlink_NORAD_IDs)

    files = [file for file in os.listdir(TLE_analysis_path) if file.endswith('.csv')]
    for file in files:
        df = process_TLE_analysis_file(file, TLE_analysis_path, oneweb_NORAD_IDs, starlink_NORAD_IDs)
        norad_id = file[:-4]
        if norad_id in oneweb_NORAD_IDs:
            oneweb_dfs.append(df)
        elif norad_id in starlink_NORAD_IDs:
            starlink_dfs.append(df)

    return oneweb_dfs, starlink_dfs

def calculate_stats(df_list):
    """Calculate the statistics for each dataframe in the list"""
    # concatenate all dataframes in the list
    df_all = pd.concat(df_list, ignore_index=True)

    # calculate statistics for each launch_no
    stats = df_all.groupby(['constellation', 'launch_no']).agg(
        mean_h_diff=('h_diffs', 'mean'),
        stddev_h_diff=('h_diffs', 'std'),
        min_h_diff=('h_diffs', 'min'),
        max_h_diff=('h_diffs', 'max'),
        mean_c_diff=('c_diffs', 'mean'),
        stddev_c_diff=('c_diffs', 'std'),
        min_c_diff=('c_diffs', 'min'),
        max_c_diff=('c_diffs', 'max'),
        mean_l_diff=('l_diffs', 'mean'),
        stddev_l_diff=('l_diffs', 'std'),
        min_l_diff=('l_diffs', 'min'),
        max_l_diff=('l_diffs', 'max'),
        mean_cart_pos_diff=('cart_pos_diffs', 'mean'),
        stddev_cart_pos_diff=('cart_pos_diffs', 'std'),
        min_cart_pos_diff=('cart_pos_diffs', 'min'),
        max_cart_pos_diff=('cart_pos_diffs', 'max')
    ).reset_index()

    print("stats: ", stats)

    return stats

def launch_specific_stats(list_of_dfs, export=True):
    """Take a lists of dataframes and return launch-specific analysis of positional differences 

    Returns:
        pd.DataFrame: Pandas dataframe containing summary statistics for each launch and each constellations
    """
    launch_summary_stats = calculate_stats(list_of_dfs)  
    print("launch summary stats: ", launch_summary_stats)

    # convert to categorical variable and int
    launch_summary_stats['launch_no'] = launch_summary_stats['launch_no'].astype('category')
    launch_summary_stats['launch_no'] = launch_summary_stats['launch_no'].astype('int')

    if export:
        launch_summary_stats.to_csv("output/launch_specific/launch_summary_stats.csv", index=False)
    return launch_summary_stats

def compute_fft(df: pd.DataFrame, diff_type: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the Fast Fourier Transform and related values on given data.

    Args:
        df (pd.DataFrame): The data frame containing the data on which FFT is to be computed.
        diff_type (str): The type of difference for which FFT is to be computed.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: The frequencies, power spectral density, 
        and an array indicating where the frequencies are greater than 0.
    """
    date = df['UTC']
    df = df.set_index('UTC')
    diff = df[diff_type]
    N = len(diff)

    diff_fft = sp.fftpack.fft(diff.values)
    diff_psd = np.abs(diff_fft)**2
    fftfreqs = sp.fftpack.fftfreq(len(diff_psd), 1./(96)) #there are 96 15-minute intervals in a day so this gives us the frequencies in days^-1
    im = fftfreqs >0

    return fftfreqs[im], 10 * np.log10(diff_psd[im])

def find_files_sup_gp_op(folder_path = 'external/ephem_TLE_compare'):
    """
    This function returns a list of file paths for files in the given directory.
    
    Args:
        folder_path (str): Path to the directory.
    
    Returns:
        list: List of file paths.
    """
    sup_list, gp_list, ephem_list = [], [], []
    for sc_folder in os.listdir(folder_path):
        for file in os.listdir(os.path.join(folder_path, sc_folder)):
            if file.startswith('sup'):
                sup_list.append(os.path.join(folder_path, sc_folder, file))
            elif file.startswith('gp'):
                gp_list.append(os.path.join(folder_path, sc_folder, file))
            elif file.startswith('MEME'):
                ephem_list.append(os.path.join(folder_path, sc_folder, file))
    return sup_list, gp_list, ephem_list

def sup_gp_op_benchmark():

    # 57.3494% age difference between sup and gp TLEs against operator ephem

    sup_list, gp_list, ephem_list = find_files_sup_gp_op()
        
#     # now going through each spacecraft
    all_triple_ephems = [] # list of the combined ephemeris data for each spacecraft
    all_sup_ages = [] # list of the ages of the suplemental TLEs
    all_gp_ages = [] # list of the ages of the gp TLEs
    all_sup_tle_epochs = [] # list of the epochs of the suplemental TLEs
    all_gp_tle_epochs = [] # list of the epochs of the gp TLEs

    for i in range(len(sup_list)):
        #read the first 5 lines of the ephem file
        sup_path = sup_list[i]
        gp_path = gp_list[i]
        ephem_path = ephem_list[i]

        # read the first 5 lines of the operator ephem file
        ephem_start_jd_dt_obj, ephem_end_jd_dt_obj, ephem_step_size = read_spacex_ephemeris(ephem_path)

        gp_TLE_list = read_TLEs(gp_path)
        sup_TLE_list = read_TLEs(sup_path)

        sup_tle_epochs = [TLE_time(TLE) for TLE in sup_TLE_list]
        all_sup_tle_epochs.append(sup_tle_epochs)

        gp_tle_epochs = [TLE_time(TLE) for TLE in gp_TLE_list]
        all_gp_tle_epochs.append(gp_tle_epochs)

        # compare the start time of the three data sources
        gp_start = TLE_time(gp_TLE_list[0])
        sup_start = TLE_time(sup_TLE_list[0])
        # if any of the TLE (sup or gp) start times are after the ephemeris start time, then return an error and end the loop
        if gp_start > ephem_start_jd_dt_obj or sup_start > ephem_start_jd_dt_obj:
            print('The start time of the TLE is after the start time of the ephemeris. Please provide a TLE that starts before the ephemeris start time.')
            break
            # this is because i cannot interpolate (using a propagator) the ephemeris, but i can interpolate the TLEs
            # so I interpolate the TLEs to the ephemeris start time
        # extract the start time of the ephemeris and set it as the start time (compare_start) for the combine_TLE2eph function
        else:
            compare_start = ephem_start_jd_dt_obj
        # now compare the end date of all data sources (ephemeris, sup, gp)
        gp_end = TLE_time(gp_TLE_list[-1])
        sup_end = TLE_time(sup_TLE_list[-1])

        # now make ephemerides using the GP and SUP TLEs
        sup_eph, sup_ages = combine_TLE2eph(TLE_list=sup_TLE_list, jd_start=ephem_start_jd_dt_obj, jd_stop=ephem_end_jd_dt_obj, dt=ephem_step_size)
        gp_eph, gp_ages = combine_TLE2eph(TLE_list=gp_TLE_list, jd_start=ephem_start_jd_dt_obj, jd_stop=ephem_end_jd_dt_obj, dt=ephem_step_size)
        
        # save the ages of the TLEs
        all_sup_ages.append(sup_ages)
        all_gp_ages.append(gp_ages)

        # now make dataframe with a row for each time step and a columns for jd, x, y, z, u,v,w
        sup_df = pd.DataFrame(columns = ['sup_jd', 'sup_x', 'sup_y', 'sup_z', 'sup_u', 'sup_v', 'sup_w'])
        gp_df = pd.DataFrame(columns = ['gp_jd', 'gp_x', 'gp_y', 'gp_z', 'gp_u', 'gp_v', 'gp_w'])

        # go through each row in sup_df and apply TEME_to_MEME to the x,y,z,u,v,w columns. Replace them inplace
        for i in range(len(sup_eph)):
            meme_x, meme_y, meme_z, meme_u, meme_v, meme_w = TEME_to_MEME(x = sup_eph[i][1][0], y = sup_eph[i][1][1], z = sup_eph[i][1][2], u = sup_eph[i][2][0], v = sup_eph[i][2][1], w = sup_eph[i][2][2], jd_time = sup_eph[i][0])
            sup_df.loc[i] = [sup_eph[i][0], meme_x, meme_y, meme_z, meme_u, meme_v, meme_w]

        for i in range(len(gp_eph)):
            meme_x, meme_y, meme_z, meme_u, meme_v, meme_w = TEME_to_MEME(x = gp_eph[i][1][0], y = gp_eph[i][1][1], z = gp_eph[i][1][2], u = gp_eph[i][2][0], v = gp_eph[i][2][1], w = gp_eph[i][2][2], jd_time = gp_eph[i][0])
            gp_df.loc[i] = [gp_eph[i][0], meme_x, meme_y, meme_z, meme_u, meme_v, meme_w]

        #merge the two dataframes
        sup_n_gp_df = pd.merge(sup_df, gp_df, left_on = 'sup_jd', right_on = 'gp_jd')

        # read in the text file 
        spacex_ephem_df = spacex_ephem_to_dataframe(ephem_path)

        # merge the two dataframes (on index since the time stamps are off by around 8 seconds after 24 hours)
        triple_ephem_df = pd.merge(sup_n_gp_df, spacex_ephem_df, left_index = True, right_index = True)

        # now calculate the H, C, L, 3D diffs for each time step between the SUP/GP and SpaceX ephemerides
        prefs = ['gp_', 'sup_']
        for pref in prefs:
            H_diffs = []
            C_diffs = []
            L_diffs = []
            cart_pos_diffs = []
            for i in range(len(triple_ephem_df)):
                tle_r = np.array([triple_ephem_df[pref + 'x'][i], triple_ephem_df[pref + 'y'][i], triple_ephem_df[pref + 'z'][i]])
                ephem_r = np.array([triple_ephem_df['x'][i], triple_ephem_df['y'][i], triple_ephem_df['z'][i]])

                tle_v = np.array([triple_ephem_df[pref + 'u'][i], triple_ephem_df[pref + 'v'][i], triple_ephem_df[pref + 'w'][i]])
                ephem_v = np.array([triple_ephem_df['u'][i], triple_ephem_df['v'][i], triple_ephem_df['w'][i]])

                unit_radial = ephem_r/np.linalg.norm(ephem_r)
                unit_cross_track = np.array(np.cross(ephem_r, ephem_v)/np.linalg.norm(np.cross(ephem_r, ephem_v)))
                unit_along_track = np.cross(unit_radial, unit_cross_track)

                unit_vectors = np.array([unit_radial, unit_cross_track, unit_along_track])

                r_diff = tle_r - ephem_r

                r_diff_HCL = np.matmul(unit_vectors, r_diff)

                h_diff = r_diff_HCL[0]
                c_diff = r_diff_HCL[1]
                l_diff = r_diff_HCL[2]

                cart_pos_diff = np.linalg.norm(tle_r - ephem_r)

                H_diffs.append(h_diff)
                C_diffs.append(c_diff)
                L_diffs.append(l_diff)
                cart_pos_diffs.append(cart_pos_diff)

            triple_ephem_df[pref + 'h_diff'] = H_diffs
            triple_ephem_df[pref + 'c_diff'] = C_diffs
            triple_ephem_df[pref + 'l_diff'] = L_diffs
            triple_ephem_df[pref + 'cart_pos_diff'] = cart_pos_diffs
            # convert 'jd_time' to mjd_time for plotting
            triple_ephem_df['mjd_time'] = triple_ephem_df['jd_time'] - 2400000.5

        all_triple_ephems.append(triple_ephem_df)
    return all_triple_ephems, all_sup_tle_epochs, all_gp_tle_epochs, gp_list

def get_NORADS_from_JSON(selected_satellites='external/selected_satellites.json'):

    # Open and load the JSON file as a Python dictionary
    with open(selected_satellites, 'r') as f:
        data = json.load(f)
    
    # Create an empty list to store the ID numbers
    ids_list = []
    
    # Iterate over the constellations
    for constellation in data.values():
        # Iterate over the levels of each constellation
        for level in constellation.values():
            # Append the ID numbers to the list
            ids_list.append(level)
    
    return ids_list

def TLE_arglat_dict(selected_satellites: str='external/selected_satellites.json', tle_folder: str='external/NORAD_TLEs/') -> List[Dict[int, List[float]]]:
    """
    Generate a list of dictionaries containing argument of latitude values for each NORAD ID in each satellite constellation.

    This function reads NORAD Two-Line Element (TLE) data for selected satellites from JSON and text files, parses and
    converts the TLE data into Keplerian elements, calculates the argument of latitude in degrees for each TLE, and 
    organizes this data into dictionaries.

    Parameters
    ----------
    selected_satellites : str, optional
        Path to the JSON file containing the selected satellites' NORAD IDs. 
        Each entry in the JSON file should correspond to a constellation. Defaults to 'external/selected_satellites.json'.
    tle_folder : str, optional
        Path to the folder containing the TLE text files for each NORAD ID. Defaults to 'external/NORAD_TLEs/'.

    Returns
    -------
    const_TLE_arglat_ls : list of dict of {int: list of float}
        List of dictionaries, where each dictionary corresponds to a constellation. 
        Each dictionary's keys are NORAD IDs, and the values are lists of argument of latitude values for each TLE.

    Notes
    -----
    This function relies on the helper functions `read_TLEs`, `twoLE_parse`, and `tle_convert`,
    which are responsible for reading the NORAD IDs from the JSON file, reading TLEs from the text files, parsing the TLEs,
    and converting the TLEs to Keplerian elements, respectively.
    """
    with open(selected_satellites, 'r') as f:
        data = json.load(f)

    const_TLE_arglat_dict = {} #dictionary of dictionaries of TLE times per NORAD ID separated by constellation

    for constellation, launches in data.items():
        TLE_arglat_dict = {} 
        for launch, NORAD_IDs in launches.items():
            for NORAD in NORAD_IDs:
                for i in os.listdir(tle_folder):
                    #keep only the text files
                    tle_file = tle_folder + i
                    sc_TLE_arglats = [] #list of TLE times for each NORAD ID
                    # get NORAD ID by getting the 5 digits before the '.txt' in the file name
                    no_txt = i[:-4]
                    #make a variable called file_NORAD that contains just the numbers from the string (no '.txt')
                    file_NORAD = int(''.join(filter(str.isdigit, no_txt)))

                    if file_NORAD == NORAD: # if the NORAD ID is part of the constellation we are looking at
                        TLEs = read_TLEs(tle_file) #read the TLEs from the text file
                        for i in range(len(TLEs)): #for each TLE
                            tle_dict = twoLE_parse(TLEs[i]) #parse the TLE
                            kep_elem = tle_convert(tle_dict)  #convert the TLE to Keplerian elements
                            # argp is in radians, true_anomaly is in radians now
                            argument_of_latitude = (kep_elem['arg_p'] + np.deg2rad(kep_elem['true_anomaly'])) % (2*np.pi) #calculate the argument of latitude
                            if argument_of_latitude > 2 * np.pi:
                                print('arglat > 2pi:', argument_of_latitude)
                            arg_lat = (argument_of_latitude * 180/np.pi) #convert the argument of latitude to degrees
                            sc_TLE_arglats.append(arg_lat) #append to the list of TLE arglats for that NORAD ID
                        TLE_arglat_dict[NORAD] = sc_TLE_arglats #add the list of TLE arglats to the dictionary for that NORAD ID
        const_TLE_arglat_dict[constellation] = TLE_arglat_dict #add the dictionary of NORAD IDs and TLE arglats to the dictionary for each constellation
    return const_TLE_arglat_dict #return the dictionary of dictionaries for each constellation

def TLE_rate_dicts(selected_satellites: str='external/selected_satellites.json', tle_folder = 'external/NORAD_TLEs/'):
    """Given a list of lists of NORAD IDs, populate a dictionary of TLE rate (hours/TLE) for each NORAD ID in each list.

    Args:
        const_norad_list (array-like): a list of lists of NORAD IDS. For example, [[NORAD_Oneweb1, NORAD_Oneweb2, NORAD_Oneweb3], [NORAD_Starlink1, NORAD_Starlink2, NORAD_Starlink3]].
        tle_folder (str, optional): Path of folder containing text files with TLE lists. Defaults to 'external/NORAD_TLEs'.

    Returns:
        list: list of dictionaries of TLE rates (hours/TLE) per NORAD ID separated by constellation.
    """

    with open(selected_satellites) as file:
        data = json.load(file)
    
    const_norad_list = []
    for constellation in data.values():
        norad_ids = []
        for norad_list in constellation.values():
            norad_ids.extend(norad_list)
        const_norad_list.append(norad_ids)

    const_tle_rates = [] # list of dictionaries of TLE rates per NORAD ID separated by constellation
    
    for const in const_norad_list:
        TLE_tdiff_dict = {}
        
        for i in os.listdir(tle_folder):
            if i[-4:] == '.txt':
                tle_file = tle_folder + i
                no_txt = i[:-4]
                NORAD = int(no_txt[-5:])
                
                if NORAD in const:
                    TLEs = read_TLEs(tle_file)
                    sc_TLE_times = [TLE_time(TLEs[i]) for i in range(len(TLEs))]  # list of TLE times for each NORAD ID
                    
                    # Calculate time differences between consecutive TLE times
                    TLE_tdiff = [(sc_TLE_times[i+1] - sc_TLE_times[i])*24 for i in range(len(sc_TLE_times) - 1)]
                    TLE_tdiff_dict[NORAD] = TLE_tdiff
        
        # Sort the dictionary by ascending NORAD ID
        TLE_tdiff_dict = dict(sorted(TLE_tdiff_dict.items()))
        const_tle_rates.append(TLE_tdiff_dict)
        
    return const_tle_rates

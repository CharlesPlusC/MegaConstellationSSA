"""
Functions to analyse the TLEs and ephemerides
"""
import os
import re
import pandas as pd
import numpy as np
import datetime

#local imports
from .tletools import read_TLEs, TLE_time, sgp4_prop_TLE, combine_TLE2eph
from .conversions import jd_to_utc, kep2car, utc_jd_date, midnight_jd_date, HCL_diff, dist_3d, alt_series

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
    
    if update == False:
        #make an empty array of the same shape as master_TLE called master_ephemeris
        master_ephemeris = np.empty_like(master_TLE, dtype=object)
        
        ## One loop that just propagates all the TLEs within the given time period
        for sc in range(0, len(master_TLE), 1): # for each spacecraft's list of TLEs
            SC_TLEs = master_TLE[sc] # list of TLEs for a single spacecraft
            
            for tle in range(0, len(SC_TLEs), 1): # for each TLE in the list of TLEs for a single spacecraft
                single_TLE = SC_TLEs[tle]
            prop = sgp4_prop_TLE(single_TLE, jd_start, jd_stop, dt) # propagate the TLE over the selected time period 
            master_ephemeris[sc][tle] = prop #fill in the relevant spot in the master ephemeris array
        
    elif update == True:
        ## Another loop that combines the ephemerides for the most recent TLEs into a single ephemeris for each spacecraft
        master_ephemeris = np.zeros(len(master_TLE), dtype=object)
        orbit_ages = np.zeros(len(master_TLE), dtype=object)
    
        for sc in range(0, len(master_TLE), 1): # for each spacecraft's list of TLEs
            SC_TLEs = master_TLE[sc] # list of TLEs for a single spacecraft
            SC_eph, SC_orbit_ages = combine_TLE2eph(SC_TLEs, jd_start, jd_stop) # combine the TLEs into a single ephemeris over the determined time period
            master_ephemeris[sc] = SC_eph
            orbit_ages[sc] = SC_orbit_ages

    return master_ephemeris, orbit_ages


def TLE_pair_analyse(pair_TLE_list, plot=False, savepath = '/home/charlesc/Documents/GitHub/Astrodynamics/source/propagation/results/images/TLE_compare/'):
    
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

    #Extract the NORAD ID of both the TLEs by parsing the first line of each of them
    norad1 = tlelist1[0][2:7]
    norad2 = tlelist2[0][2:7]

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

    print('Start date: ', start_year, start_month, start_day)
    print('End date: ', end_year, end_month, end_day)

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

        cart_rmse = rmse(cart_pos_diffs)

    # Project into HCL 
    h_diffs, c_diffs, l_diffs = HCL_diff(pos_vels[0], pos_vels[1], plot=False)

    #altitude time series for both 
    eph_alts = []
    for i in range(0, len(master_ephs), 1):
        alts = alt_series(master_ephs[i])
        eph_alts.append(alts)
        
        return master_ephs,eph_alts, h_diffs, c_diffs, l_diffs, cart_pos_diffs, times, orbit_ages

def NORAD_vs_SUP_TLE_analysis(NORADS = [], analysis_output_path = 'output/TLE_analysis'):

    """Analyzes the quality of NORAD and Operator based orbits (propagated using SGP4).
       Will search each folder for mathching NORAD IDs and then compare the TLEs by propagating them using SGP4, and updating them as soon as a new TLE is available.
       Every 15 minutes, the height cross track and along track differences are recorded. The ephemerides are also recorded.

    Args:
        TLE_source1_path (str): path to a folder containing TLEs from the first source
        TLE_source2_path (str): path to a folder containing TLEs from the second source
        analysis_output_path (str): path to a folder where the analysis will be saved. One CSV file will be created for each NORAD ID.
     """

    NORAD_TLE_folder = 'external/NORAD_TLEs/'
    SUP_TLE_folder = 'external/SUP_TLEs/'

    NORAD_TLE_files = os.listdir(NORAD_TLE_folder)
    SUP_TLE_files = os.listdir(SUP_TLE_folder)

    if NORADS == []:
        raise ValueError('Please specify a list of NORAD IDs to analyze')
    else:
        #if a specific NORAD ID is specified, then find the norad TLE files that correspond to those NORAD IDs
        NORAD_TLE_files = [file for file in NORAD_TLE_files if file.startswith(tuple(NORADS))]
        # find the files in SUP_TLE_files that correspond to the NORAD IDs with a 'SUP' prefix
        SUP_TLE_files = [file for file in SUP_TLE_files if file.startswith(tuple(['sup_' + NORAD for NORAD in NORADS]))]
        #find the files with the same NORAD ID

        for NORAD_file in NORAD_TLE_files:
            for SUP_file in SUP_TLE_files:
                #extract the numerical values from NORAD_file and SUP_file
                NORAD_id = re.findall(r'\d+', NORAD_file)
                SUP_id = re.findall(r'\d+', SUP_file) 
                NORAD_NORAD = NORAD_id
                SUP_NORAD = SUP_id
                if NORAD_NORAD == SUP_NORAD:
                    print('Match found for NORAD ID:', NORAD_NORAD)
                    filename = str(NORAD_NORAD[0]) + '.csv'
                    total_out_path = analysis_output_path + filename
                    #set the paths to the files
                    NORAD_TLE_path = NORAD_TLE_folder + NORAD_file
                    SUP_TLE_path = SUP_TLE_folder + SUP_file
                    #Read the TLEs
                    sup_read = read_TLEs(SUP_TLE_path)
                    print('SUP TLEs read')
                    NORAD_read = read_TLEs(NORAD_TLE_path)
                    print('NORAD TLEs read')

                    #combine the TLEs into a list 
                    sup_NORAD_pair = [sup_read, NORAD_read]
                    
                    #Run the comparison function
                    master_ephs, eph_alts, h_diffs, c_diffs, l_diffs, cart_pos_diffs, times, orbit_ages = TLE_pair_analyse(sup_NORAD_pair)

                    # Make pandas dataframe and save into csv
                    df = pd.DataFrame({'h_diffs': h_diffs, 'c_diffs': c_diffs, 'l_diffs': l_diffs, 'cart_pos_diffs': cart_pos_diffs, 'times': times, 'eph_alts_sup': eph_alts[0], 'eph_alts_norad': eph_alts[1], 'master_ephs_sup': master_ephs[0], 'master_ephs_norad': master_ephs[1], 'orbit_ages_sup': orbit_ages[0], 'orbit_ages_norad': orbit_ages[1]})

                    df.to_csv(total_out_path)
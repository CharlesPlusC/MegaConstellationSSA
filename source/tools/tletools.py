""" 
Functions to download, propagate and manipulate TLE data. 
"""

#imports
import os
import configparser
import json
import time
import requests
import numpy as np
from sgp4.api import Satrec
from datetime import datetime, timedelta

#local imports
from .conversions import kep2car

class MyError(Exception):
    def __init___(self, args):
        Exception.__init__(
            self, "my exception was raised with arguments {0}".format(args)
        )
        self.args = args

def load_satellite_lists(file_path="external/selected_satellites.json"):
    with open(file_path, "r") as f:
        satellite_lists = json.load(f)
    return satellite_lists

def SpaceTrack_authenticate():
    """Authenticate with SpaceTrack using the credentials stored in the config file. 

    Returns:
        str: SpaceTrack username
        str: SpaceTrack password

    """
    #if you set environment variables these can be fetched automatically here
    password = os.environ.get("SLTRACK_PWD")
    username = os.environ.get("SLTRACK_USR")
    # if they are not found, then it will look in SLTrack.ini
    if not all((username, password)):

        config = configparser.ConfigParser()
        config.read("SLTrack.ini")
        username = config.get("configuration", "username")
        password = config.get("configuration", "password")

    assert ( username is not None and password is not None), "Please specify a username and password for SpaceTrack"
    
    return username, password

def get_satellite_ids(session, uriBase, requestCmdAction, requestFindStarlinks):
    resp = session.get(uriBase + requestCmdAction + requestFindStarlinks)
    if resp.status_code != 200:
        print(resp)
        raise Exception("GET fail on request for satellites")

    retData = json.loads(resp.text)
    return [int(e['NORAD_CAT_ID']) for e in retData]

def wait_for_next_request(lastRequest):
    if time.time() - lastRequest < 13:
        time.sleep(13 - (time.time() - lastRequest))

def get_satellite_data(session, uriBase, requestCmdAction, requestOMMStarlink1, requestOMMStarlink2, s):
    wait_for_next_request(time.time())
    resp = session.get(uriBase + requestCmdAction + requestOMMStarlink1 + str(s) + requestOMMStarlink2)
    if resp.status_code != 200:
        print(resp)
        raise Exception("GET fail on request for satellite " + s)

    return resp.text

def process_satellite_data(output, folder_path, s):
    # Process the TLE data and write it to the file that is named after its NORAD ID
    with open(folder_path + str(s) + '.txt', "w") as f:
        all_lines = output.splitlines()
        stringofdicts = all_lines[0]
        dicts = stringofdicts.split('},')
        dicts = [x[2:] for x in dicts]
        dicts = [x.replace('"', '') for x in dicts]
        for i in range(len(dicts)):
            key_value_pairs = dicts[i].split(',')
            key_value_pairs = [x.split(':') for x in key_value_pairs]
            spacetrack_dict = dict(key_value_pairs)

            # this will write the pulled TLE to the text file
            f.write(spacetrack_dict['TLE_LINE1'] + '\n' + spacetrack_dict['TLE_LINE2'] + '\n')
        f.close()# Process the TLE data and write it to the file

def tle_convert(tle_dict, display=False):
    """
    Converts a TLE dictionary into the corresponding keplerian elements
    
    Args:
        tle_dict (dict): dictionary of TLE data as provided by the tle_parse function

    Returns:
        keplerian_dict(dict): dictionary containing Keplerian elements
    """

    # Standard gravitational parameter for the Earth
    GM = 398600.4415 * (1e3)**3 # m^3/s^2

    # Convert RAAN from degrees to radians
    RAAN = np.radians(float(tle_dict['right ascension of the ascending node']))
    
    # Convert argument of perigee from degrees to radians
    arg_p = np.radians(float(tle_dict['argument of perigee']))
    
    # Convert mean motion from revolutions per day to radians per second
    mean_motion = float(tle_dict['mean motion']) * (2 * np.pi / 86400)
    
    # Compute the period of the orbit in seconds
    period = 2 * np.pi / mean_motion
    
    # Compute the semi-major axis
    n = mean_motion # mean motion in radians per second
    a = (GM / (n ** 2)) ** (1/3) / 1000 # in km
    
    # Convert mean anomaly from degrees to radians
    M = np.radians(float(tle_dict['mean anomaly']))
    
    # Extract eccentricity as decimal value
    e = float("0." + tle_dict['eccentricity'])
    
    # Convert inclination from degrees to radians
    inclination = np.radians(float(tle_dict['inclination']))
    
    # Initial Guess at Eccentric Anomaly
    if M < np.pi:
        E = M + (e / 2)
    else:
        E = M - (e / 2)

    # Numerical iteration for Eccentric Anomaly
    f = lambda E: E - e * np.sin(E) - M
    fp = lambda E: 1 - e * np.cos(E)
    E = np.float64(E)
    r_tol = 1e-8 # set the convergence tolerance for the iteration
    max_iter = 50 # set the maximum number of iterations allowed
    for it in range(max_iter):
        f_value = f(E)
        fp_value = fp(E)
        E_new = E - f_value / fp_value
        if np.abs(E_new - E) < r_tol:
            E = E_new
            break
        E = E_new
    else:
        raise ValueError("Eccentric anomaly did not converge")
        
    eccentric_anomaly = E

    # Compute True Anomaly
    true_anomaly = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(eccentric_anomaly / 2),
                                  np.sqrt(1 - e) * np.cos(eccentric_anomaly / 2))

    # Dictionary of Keplerian elements
    keplerian_dict = {'a': a, 'e': e, 'i': inclination, 'RAAN': RAAN, 'arg_p': arg_p, 'true_anomaly': np.degrees(true_anomaly)}
    if display == True:
        print("Keplerian Elements:")
        print("a = {:.2f} km".format(keplerian_dict['a']))
        print("e = {:.2f}".format(keplerian_dict['e']))
        print("i = {:.2f} deg".format(np.degrees(keplerian_dict['i'])))
        print("RAAN = {:.2f} deg".format(np.degrees(keplerian_dict['RAAN'])))
        print("arg_p = {:.2f} deg".format(np.degrees(keplerian_dict['arg_p'])))
        print("true_anomaly = {:.2f} deg".format(np.degrees(keplerian_dict['true_anomaly'])))

    return keplerian_dict

def download_tle_history(NORAD_ids, constellation, folder_path="external/NORAD_TLEs"):
    """
    This function takes a list of NORAD IDs and returns a dictionary of all available TLEs for each satellite. 
    The TLEs are returned as a list of strings, with each string containing the TLE for a single satellite. The dictionary keys are the NORAD IDs.
    The TLEs are returned in chronological order, with the most recent TLE first. The function samples the entire archive of available TLEs for each satellite.
    NOTE: The function will take a long time to run if there are many satellites in the list because there is a rate limit on the number of queries that can be made to the API (1 query per 13 seconds).

    Args:
        NORAD_ids (list): A list of NORAD IDs for the satellites of interest.
        constellation (str): The name of the constellation. This is used to determine the correct TLE archive to query.
        folder_path (str): The path to the folder where the TLEs will be saved. Defaults to "external/NORAD_TLEs".
    """

    available_constellations = ["starlink", "oneweb"]
    if constellation not in available_constellations:
        raise ValueError(f"Invalid constellation name. Select one of: {available_constellations}")

    if constellation == "starlink":
        const = "STARLINK"
    elif constellation == "oneweb":
        const = "ONEWEB"

    # Verify that NORAD IDs are integers
    if not all(isinstance(NORAD_id, int) for NORAD_id in NORAD_ids):
        raise ValueError("NORAD ids must all be integers")

    # Authenticate with Space-Track
    username, password = SpaceTrack_authenticate()
    site_credentials = {"identity": username, "password": password}

    # Define API endpoints
    uri_base = "https://www.space-track.org"
    login_endpoint = "/ajaxauth/login"
    query_endpoint = "/basicspacedata/query"
    find_satellites_query = f"/class/tle_latest/NORAD_CAT_ID/>40000/ORDINAL/1/OBJECT_NAME/{const}~~/format/json/orderby/NORAD_CAT_ID%20asc"
    omm_query_1 = "/class/omm/NORAD_CAT_ID/"
    omm_query_2 = "/orderby/EPOCH%20asc/format/json"

    # Use requests package to drive the RESTful session with space-track.org
    with requests.Session() as session:
        # Log in to Space-Track
        response = session.post(uri_base + login_endpoint, data=site_credentials)
        if response.status_code != 200:
            raise ValueError(f"POST failed on login: {response}")

        # Query for satellites
        response = session.get(uri_base + query_endpoint + find_satellites_query)
        if response.status_code != 200:
            raise ValueError(f"GET failed on request for satellites: {response}")

        # Parse response JSON
        satellite_data = json.loads(response.text)
        satellite_ids = [int(entry["NORAD_CAT_ID"]) for entry in satellite_data]

        # Download TLE history for each satellite
        for norad_id in NORAD_ids:
            if norad_id in satellite_ids:
                # Throttle requests to avoid rate limiting (1 query per 13 seconds)
                last_request = time.time()
                if time.time() - last_request < 13:
                    time.sleep(13 - (time.time() - last_request))
                last_request = time.time()

                # Get OMM data for the specific NORAD ID
                response = session.get(uri_base + query_endpoint + omm_query_1 + str(norad_id) + omm_query_2)
                if response.status_code != 200:
                    raise ValueError(f"GET failed on request for satellite {norad_id}: {response}")

                # Save TLE data to a text file
                output = response.text
                with open(f"{folder_path}/{norad_id}.txt", "w") as f:
                    # Split the output into lines
                    all_lines = output.splitlines()
                    #
                    list_of_dicts = all_lines[0]
                    data_list = json.loads(list_of_dicts)
                    # Write the TLEs to the file
                    for line in data_list:
                        f.write(line["TLE_LINE1"] + "\n")
                        f.write(line["TLE_LINE2"] + "\n")

                print(f"Downloaded TLE history for satellite {norad_id}")
            else:
                print(f"Satellite {norad_id} not found in specified constellation: {constellation}")

def read_TLEs(filename):
    """Read a TLE file and return a list of TLEs

    Args:
        filename (string): name of the TLE file

    Returns:
        list: list of TLEs
    """
    #open the file
    with open(filename, 'r') as f:
        #read the file
        TLEs = f.readlines()
        #split the file into a list of TLEs every 2 lines
        TLEs = [TLEs[i:i+2] for i in range(0, len(TLEs), 2)]
        #remove the new line character from the end of each TLE
        TLEs = [TLEs[i][0] + '' + TLEs[i][1].strip('\n') for i in range(0, len(TLEs), 1)]
        #drop the last two characters of each TLE (get rid of the \n)
        print("number of TLEs read: ", len(TLEs))
        #close the file
        f.close()
        #return the list of TLEs
        return TLEs
    
def TLE_time(TLE):
    """Find the time of a TLE in julian day format"""
    #find the epoch section of the TLE
    epoch = TLE[18:32]
    #convert the first two digits of the epoch to the year
    year = 2000+int(epoch[0:2])
    
    # the rest of the digits are the day of the year and fractional portion of the day
    day = float(epoch[2:])
    #convert the day of the year to a day, month, year format
    date = datetime.datetime(year, 1, 1) + datetime.timedelta(day - 1)
    #convert the date to a julian date
    jd = (date - datetime.datetime(1858, 11, 17)).total_seconds() / 86400.0 + 2400000.5
    return jd

def sgp4_prop_TLE(TLE, jd_start, jd_end, dt, alt_series = False):

    """Given a TLE, a start time, end time, and time step, propagate the TLE and return the time-series of Cartesian coordinates, and accompanying time-stamps (MJD)
        Note: Simply a wrapper for the SGP4 routine in the sgp4.api package (Brandon Rhodes)
    Args:
        TLE (string): TLE to be propagated
        jd_start (float): start time of propagation in Julian Date format
        jd_end (float): end time of propagation in Julian Date format
        dt (float): time step of propagation in seconds
        alt_series (bool, optional): If True, return the altitude series as well as the position series. Defaults to False.
        
    Returns:
    list: list of lists containing the time-series of Cartesian coordinates, and accompanying time-stamps (MJD)
    
    """

    if jd_start > jd_end:
        print('jd_start must be less than jd_end')
        return

    ephemeris = []
    
    #convert dt from seconds to julian day
    dt_jd = dt/86400

    #split at the new line
    split_tle = TLE.split('\n')
    s = split_tle[0]
    r = split_tle[1]

    fr = 0.0 # precise fraction (SGP4 docs for more info)
    
    #create a satellite object
    satellite = Satrec.twoline2rv(s, r)

    time = jd_start
    # for i in range (jd_start, jd_end, dt):
    while time < jd_end:
        # propagate the satellite to the next time step
        # Position is in idiosyncratic True Equator Mean Equinox coordinate frame used by SGP4
        # Velocity is the rate at which the position is changing, expressed in kilometers per second
        error, position, velocity = satellite.sgp4(time, fr)
        if error != 0:
            print('Satellite position could not be computed for the given date')
            break
        else:
            ephemeris.append([time,position, velocity]) #jd time, pos, vel
        time += dt_jd

    return ephemeris


def combine_TLE2eph(TLE_list,jd_start,jd_stop, dt=(15*60)):
    """ Takes a list of TLES a returns an ephemeris that updates with each new TLE. Outputs a position and velocity every 15 minutes from the hour.
    Args:
        TLE_list (list): list of TLEs (use read_TLEs function to generate this)
        jd_start (float): start time in JD
        jd_stop (float): stop time in JD
        dt (float): time step in seconds
    Returns:
        ephemeris (array): ephemeris of the satellite in ECI coordinates(time, pos, vel)
        """
    #Given a list of TLEs return an ephemeris that is updated with the TLE
    #default dt is 15 minutes
    #propagated with sgp4
    dt_jd = dt/86400
    current_jd = jd_start #this should be the midnight of the day you want to start
    #number of steps
    n_steps = int((jd_stop - jd_start)/dt_jd)
    #t since update
    orbit_ages = []
    ephemeris = []
    while current_jd < jd_stop:
        #loop through the TLEs
        for i in range(0, len(TLE_list), 1):
            #get the time of the current TLE
            TLE_jd = TLE_time(TLE_list[i])
            #get the time of the next TLE
            if i == len(TLE_list)-1:
                next_TLE_jd = TLE_time(TLE_list[0])
            else:
                next_TLE_jd = TLE_time(TLE_list[i+1])
            #if the current jd is between the current TLE and the next TLE, use the current TLE
            if current_jd > TLE_jd and current_jd < next_TLE_jd:
                #get the ephemeris of the satellite
                eph = sgp4_prop_TLE(TLE_list[i], current_jd, (current_jd+dt_jd), dt=dt) #prop for one step 
                #append the contents of eph to the ephemeris
                ephemeris.extend(eph)
                #increment the current jd
                current_jd += dt_jd
                #increment the time since update
                jd_orbit_age = current_jd - TLE_jd
                # convert from julian days to hours
                hours_orbit_age = jd_orbit_age*24
                orbit_ages.append(hours_orbit_age)
            elif current_jd > jd_stop:
                print("prop time is greater than stop time. Stopping propagation.")
            break  

    # chop the ephemeris to be the correct number of steps using the n_steps variable (stops lists being 1 too long due to Python indexing)
    ephemeris = ephemeris[0:n_steps]
    orbit_ages = orbit_ages[0:n_steps]
    return ephemeris, orbit_ages

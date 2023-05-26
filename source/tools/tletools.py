""" Functions to load, propagate and manipulate TLE data. """
#imports
import os
import configparser
import json
import time
import requests
import numpy as np
from sgp4.api import Satrec
import datetime
from astropy.time import Time
import pandas as pd
from typing import Dict, List, Union, Tuple

#local imports
from .conversions import parse_spacex_datetime_stamps, yyyy_mm_dd_hh_mm_ss_to_jd

class MyError(Exception):
    def __init__(self, args: str):
        """
        Custom exception for satellite tracking errors.

        Parameters
        ----------
        args : str
            Arguments or error message to pass to exception.
        """
        super().__init__(
            f"my exception was raised with arguments {args}"
        )
        self.args = args

def load_satellite_lists(file_path: str = "external/selected_satellites.json") -> dict:
    """
    Load the list of satellites from a json file.

    Parameters
    ----------
    file_path : str, optional
        The path to the file containing the list of satellites, by default "external/selected_satellites.json"

    Returns
    -------
    dict
        The list of satellites.
    """
    with open(file_path, "r") as f:
        satellite_lists = json.load(f)
    return satellite_lists

def SpaceTrack_authenticate() -> Tuple[str, str]:
    """
    Authenticate with SpaceTrack using the credentials stored in the config file. 

    Returns
    -------
    Tuple[str, str]
        SpaceTrack username and password
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

def get_satellite_ids(session: requests.Session, uriBase: str, requestCmdAction: str, requestFindStarlinks: str) -> List[int]:
    """
    Get satellite ids from session.

    Parameters
    ----------
    session : requests.Session
        Session of request.
    uriBase : str
        Base uri.
    requestCmdAction : str
        Command action for the request.
    requestFindStarlinks : str
        Find Starlinks request.

    Returns
    -------
    List[int]
        List of satellite ids.

    Raises
    ------
    Exception
        If the GET request fails.
    """
    resp = session.get(uriBase + requestCmdAction + requestFindStarlinks)
    if resp.status_code != 200:
        print(resp)
        raise Exception("GET fail on request for satellites")

    retData = json.loads(resp.text)
    return [int(e['NORAD_CAT_ID']) for e in retData]

def wait_for_next_request(lastRequest: float) -> None:
    """
    Wait for next request.

    Parameters
    ----------
    lastRequest : float
        Timestamp of the last request.

    """
    if time.time() - lastRequest < 13:
        time.sleep(13 - (time.time() - lastRequest))

def get_satellite_data(session: requests.Session, uriBase: str, requestCmdAction: str, requestOMMStarlink1: str, requestOMMStarlink2: str, s: str) -> str:
    """
    Makes a GET request for a specific satellite's data.

    Parameters
    ----------
    session : requests.Session
        The active Session object.
    uriBase : str
        The base URL for the SpaceTrack API.
    requestCmdAction : str
        The command action string for the SpaceTrack API request.
    requestOMMStarlink1 : str
        The first part of the specific request string for the satellite data.
    requestOMMStarlink2 : str
        The second part of the specific request string for the satellite data.
    s : str
        The satellite identifier.

    Returns
    -------
    str
        The response text from the GET request.
    """
    wait_for_next_request(time.time())
    resp = session.get(uriBase + requestCmdAction + requestOMMStarlink1 + str(s) + requestOMMStarlink2)
    if resp.status_code != 200:
        print(resp)
        raise Exception("GET fail on request for satellite " + s)

    return resp.text

def NORAD_list_update(constellation, out_path = 'external/Constellation_NORAD_IDs/'):
    """
    From spacetrack API, get the list of all available NORAD IDs for a given constellation. Add them to a text file

    Parameters
    ----------
    constellation : str
        Constellation name. Select from: 'starlink', 'oneweb'.
    out_path : str, optional
        Path to output textfile to. Defaults to None. If None, will use default Constellation_NORAD_IDs folder.

    Raises
    ------
    ValueError
        Invalid constellation name. Select one of: ['starlink', 'oneweb']

    Returns
    -------
    list
        List of the NORAD IDs that were returned from query 

    Notes
    -----
    See https://www.space-track.org/documentation for details on REST queries
    """

    available_constellations = ['starlink','oneweb']
    if constellation not in available_constellations:
        raise ValueError("Invalid constellation name. Select one of: %s" % available_constellations)

    #navigate to the folder
    if constellation == 'oneweb':
        out_path = '/home/charlesc/Documents/GitHub/Astrodynamics/source/propagation/results/data/Constellation_NORAD_IDs/oneweb_NORAD_IDs.txt'
    elif constellation == 'starlink':
        out_path = '/home/charlesc/Documents/GitHub/Astrodynamics/source/propagation/results/data/Constellation_NORAD_IDs/starlink_NORAD_IDs.txt'
    print("NORAD list output path: ", out_path)

    uriBase = "https://www.space-track.org"
    requestLogin = "/ajaxauth/login"
    requestCmdAction = "/basicspacedata/query"
    if constellation == 'oneweb':
        requestFindOWs = "/class/tle_latest/NORAD_CAT_ID/>40000/ORDINAL/1/OBJECT_NAME/ONEWEB~~/format/tle/orderby/NORAD_CAT_ID%20asc"
    elif constellation == 'starlink':
        requestFindOWs = "/class/tle_latest/NORAD_CAT_ID/>40000/ORDINAL/1/OBJECT_NAME/Starlink~~/format/tle/orderby/NORAD_CAT_ID%20asc"

    # Find credentials for the SpaceTrack API
    
    username, password = SpaceTrack_authenticate()
    siteCred = {'identity': username, 'password': password}

    with open(out_path, "w") as f:
        # use requests package to drive the RESTful
        # session with space-track.org

        with requests.Session() as session:
            # run the session in a with block to
            # force session to close if we exit

            # need to log in first. note that we get a 200
            # to say the web site got the data, not that we are logged in
            resp = session.post(uriBase + requestLogin, data=siteCred)
            if resp.status_code != 200:
                raise MyError(resp, "POST fail on login")

            # this query picks up all OneWeb satellites from the catalog.
            # Note - a 401 failure shows you have bad credentials
            resp = session.get(uriBase + requestCmdAction + requestFindOWs)
            if resp.status_code != 200:
                raise MyError(
                    resp, "GET fail on request for satellites"
                )

            output = resp.text

        NORAD_ids = []
        f.write("last updated: " + str(datetime.now()) + "\n")
        
        # split the output into lines
        all_lines = output.splitlines()
        # put every two lines into a list
        line_pairs = [
            all_lines[i:i + 2] for i in range(0, len(all_lines), 2)
        ]

        for tle in range(0, len(line_pairs), 1):
            line_one, line_two = line_pairs[tle][0], line_pairs[tle][1]

            tle_dict = {}

            # Parse the first line
            tle_dict["line number"] = line_one[0]
            tle_dict["satellite catalog number"] = line_one[2:7]
            f.write(tle_dict["satellite catalog number"] + "\n")
            NORAD_ids.append(tle_dict["satellite catalog number"])
    return NORAD_ids

def process_satellite_data(output: str, folder_path: str, s: str) -> None:
    """
    Processes satellite TLE data and writes it to a file.

    Parameters
    ----------
    output : str
        The response text from the GET request containing TLE data.
    folder_path : str
        The path to the folder where the TLE data will be written.
    s : str
        The satellite identifier.
    """
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

def tle_convert(tle_dict: dict) -> dict:
    """
    Converts a TLE dictionary into the corresponding Keplerian elements.

    Parameters
    ----------
    tle_dict : dict
        Dictionary of TLE data.

    Returns
    -------
    dict
        Dictionary containing Keplerian elements.
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
    return keplerian_dict

def twoLE_parse(tle_2le: str) -> dict:
    """
    Parses a 2-line element set (2LE) string and returns the data in a dictionary.

    Parameters
    ----------
    tle_2le : str
        The 2-line element set string to be parsed.

    Returns
    -------
    dict
        Dictionary of the data contained in the TLE string.
    """

    # This function takes a TLE string and returns a dictionary of the TLE data
    tle_lines = tle_2le.split('\n')
    tle_dict = {}
    line_one, line_two = tle_lines[0],tle_lines[1]
    
    #Parse the first line
    tle_dict['line number'] = line_one[0]
    tle_dict['satellite catalog number'] = line_one[2:7]
    tle_dict['classification'] = line_one[7]
    tle_dict['International Designator(launch year)'] = line_one[9:11] 
    tle_dict['International Designator (launch num)'] = line_one[11:14]
    tle_dict['International Designator (piece of launch)'] = line_one[14:17]
    tle_dict['epoch year'] = line_one[18:20]
    tle_dict['epoch day'] = line_one[20:32]
    tle_dict['first time derivative of mean motion(ballisitc coefficient)'] = line_one[33:43]
    tle_dict['second time derivative of mean motion(delta-dot)'] = line_one[44:52]
    tle_dict['bstar drag term'] = line_one[53:61]
    tle_dict['ephemeris type'] = line_one[62]
    tle_dict['element number'] = line_one[63:68]
    tle_dict['checksum'] = line_one[68:69]

    #Parse the second line (ignore the line number, satellite catalog number, and checksum)
    tle_dict['inclination'] = line_two[8:16]
    tle_dict['right ascension of the ascending node'] = line_two[17:25]
    tle_dict['eccentricity'] = line_two[26:33]
    tle_dict['argument of perigee'] = line_two[34:42]
    tle_dict['mean anomaly'] = line_two[43:51]
    tle_dict['mean motion'] = line_two[52:63]
    tle_dict['revolution number at epoch'] = line_two[63:68]

    return tle_dict

def download_tle_history(NORAD_ids: List[int], constellation: str, folder_path: str = "external/NORAD_TLEs") -> None:
    """
    Downloads TLE history for a given list of NORAD IDs.
    Note: This function takes a long time to run for large satellite lists due to API rate limits.

    Parameters
    ----------
    NORAD_ids : list of int
        List of NORAD IDs for the satellites of interest.
    constellation : str
        Name of the constellation. Used to determine the correct TLE archive to query.
    folder_path : str, optional
        Path to the folder where the TLEs will be saved. Defaults to "external/NORAD_TLEs".
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

def read_TLEs(filename: str) -> List[str]:
    """
    Read a TLE file and return a list of TLEs.

    Parameters
    ----------
    filename : str
        Name of the TLE file.

    Returns
    -------
    list
        List of TLEs.
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
    
def TLE_time(TLE: str) -> float:
    """
    Find the time of a TLE in Julian Day format.

    Parameters
    ----------
    TLE : str
        The TLE string.

    Returns
    -------
    float
        Time in Julian Day format.
    """
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

def sgp4_prop_TLE(TLE: str, jd_start: float, jd_end: float, dt: float, alt_series: bool = False) -> List[List[Any]]:
    """
    Given a TLE, a start time, end time, and time step, propagate the TLE and return the time-series of Cartesian coordinates,
    and accompanying time-stamps (MJD).
    Note: Simply a wrapper for the SGP4 routine in the sgp4.api package (Brandon Rhodes).

    Parameters
    ----------
    TLE : str
        TLE to be propagated.
    jd_start : float
        Start time of propagation in Julian Date format.
    jd_end : float
        End time of propagation in Julian Date format.
    dt : float
        Time step of propagation in seconds.
    alt_series : bool, optional
        If True, return the altitude series as well as the position series. Defaults to False.

    Returns
    -------
    list
        List of lists containing the time-series of Cartesian coordinates, and accompanying time-stamps (MJD).
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

def combine_TLE2eph(TLE_list: List[str], jd_start: float, jd_stop: float, dt: float=(15 * 60)) -> Tuple[List[Any], List[Any]]:
    """
    Takes a list of TLEs and returns an ephemeris that updates with each new TLE. Outputs a position and velocity every 15 minutes from the hour.

    Parameters
    ----------
    TLE_list : list
        List of TLEs (use read_TLEs function to generate this).
    jd_start : float
        Start time in JD.
    jd_stop : float
        Stop time in JD.
    dt : float
        Time step in seconds.

    Returns
    -------
    Tuple[List[Any], List[Any]]
        Ephemeris of the satellite in ECI coordinates(time, pos, vel) and orbit ages.
    """
    dt_jd = dt / 86400
    current_jd = jd_start
    n_steps = int((jd_stop - jd_start) / dt_jd)
    ephemeris = []
    orbit_ages = []

    # Keep track of the current TLE index
    current_tle_idx = 0

    while current_jd < jd_stop:
        found_tle = False  # Flag to track if a matching TLE is found

        for i in range(current_tle_idx, len(TLE_list)):
            TLE_jd = TLE_time(TLE_list[i])
            next_TLE_jd = TLE_time(TLE_list[i + 1]) if i < len(TLE_list) - 1 else TLE_time(TLE_list[0])

            if TLE_jd < current_jd < next_TLE_jd:
                eph = sgp4_prop_TLE(TLE_list[i], current_jd, (current_jd + dt_jd), dt=dt)
                ephemeris.extend(eph)
                current_jd += dt_jd
                hours_orbit_age = (current_jd - TLE_jd) * 24
                orbit_ages.append(hours_orbit_age)
                current_tle_idx = i  # Update the TLE index
                found_tle = True
                break

        if not found_tle:
            break  # Break out of the outer loop if no matching TLE is found

    ephemeris = ephemeris[:n_steps]
    orbit_ages = orbit_ages[:n_steps]

    return ephemeris, orbit_ages

def read_spacex_ephemeris(ephem_path: str) -> Tuple[float, float, int]:
    """
    Reads a SpaceX ephemeris file and extracts start time, end time, and step size.

    Parameters
    ----------
    ephem_path : str
        Path to the ephemeris file.

    Returns
    -------
    Tuple[float, float, int]
        Start time (JD), end time (JD), step size (seconds).
    """
    # read the first 5 lines of the operator ephem file
    with open(ephem_path) as f:
        ephem_lines = f.readlines()
    ephem_lines = ephem_lines[:5]
    ephem_utc_start = str(ephem_lines[1][16:16+19]) # start time
    ephem_utc_end = str(ephem_lines[1][55:55+19]) # end time
    ephem_step_size = int(ephem_lines[1][89:89+2]) # step size
    #convert to datetime object
    ephem_utc_dt_obj_start = datetime.datetime.strptime(ephem_utc_start, '%Y-%m-%d %H:%M:%S')
    ephem_utc_dt_obj_end = datetime.datetime.strptime(ephem_utc_end, '%Y-%m-%d %H:%M:%S')
    # convert to julian date
    ephem_start_jd_dt_obj = Time(ephem_utc_dt_obj_start).jd
    ephem_end_jd_dt_obj = Time(ephem_utc_dt_obj_end).jd

    return ephem_start_jd_dt_obj, ephem_end_jd_dt_obj, ephem_step_size

def spacex_ephem_to_dataframe(ephem_path: str) -> pd.DataFrame:
    """
    Converts SpaceX ephemeris data into a pandas DataFrame. 

    Parameters
    ----------
    ephem_path : str
        Path to the ephemeris file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing parsed SpaceX ephemeris data.
    """
    # read in the text file 
    with open(ephem_path) as f:
        lines = f.readlines()
    # remove the header lines
    lines = lines[4:]
    # select every 4th line
    t_xyz_uvw = lines[::4]
    # from all the lines in t_xyz_uvw select the first float in each line and append that to a list
    t = [float(i.split()[0]) for i in t_xyz_uvw]
    x = [float(i.split()[1]) for i in t_xyz_uvw]
    y = [float(i.split()[2]) for i in t_xyz_uvw]
    z = [float(i.split()[3]) for i in t_xyz_uvw]
    u = [float(i.split()[4]) for i in t_xyz_uvw]
    v = [float(i.split()[5]) for i in t_xyz_uvw]
    w = [float(i.split()[6]) for i in t_xyz_uvw]
    
    # make all the values in the list 't' into a numpy array
    tstamp_array = np.array(t)
    # parse the timestamps into year, day of year, hour, minute, and second
    parsed_tstamps = parse_spacex_datetime_stamps(tstamp_array)
    # convert the parsed timestamps into julian dates
    jd_stamps = np.zeros(len(parsed_tstamps))
    for i in range(0, len(parsed_tstamps), 1):
        jd_stamps[i] = yyyy_mm_dd_hh_mm_ss_to_jd(int(parsed_tstamps[i][0]), int(parsed_tstamps[i][1]), int(parsed_tstamps[i][2]), int(parsed_tstamps[i][3]), int(parsed_tstamps[i][4]), int(parsed_tstamps[i][5]), int(parsed_tstamps[i][6]))

    # take t, x, y, z, u, v, w and put them into a dataframe
    spacex_ephem_df = pd.DataFrame({'jd_time':jd_stamps, 'x':x, 'y':y, 'z':z, 'u':u, 'v':v, 'w':w})
    # use the function meme_2_teme() to convert the x, y, z, u, v, w values from the MEME frame to the TEME frame
    # remove the last row from spacex_ephem_df
    spacex_ephem_df = spacex_ephem_df[:-1]

    return spacex_ephem_df
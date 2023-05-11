#imports
import os
import configparser
import json
import time
import requests
import numpy as np
from datetime import datetime, timedelta

#local imports
from ..tools.coordinates import kep2car

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

def NORAD_TLE_History(NORAD_ids, constellation, folder_path="/external/NORAD_TLEs/"):
    available_constellations = ['starlink', 'oneweb']
    if constellation not in available_constellations:
        raise ValueError("Invalid constellation name. Select one of: %s" % available_constellations)

    const = 'STARLINK' if constellation == 'starlink' else 'ONEWEB'

    uriBase = "https://www.space-track.org"
    requestLogin = "/ajaxauth/login"
    requestCmdAction = "/basicspacedata/query"
    requestFindStarlinks = "/class/tle_latest/NORAD_CAT_ID/>40000/ORDINAL/1/OBJECT_NAME/"+const+"~~/format/json/orderby/NORAD_CAT_ID%20asc"
    requestOMMStarlink1 = "/class/omm/NORAD_CAT_ID/"
    requestOMMStarlink2 = "/orderby/EPOCH%20asc/format/json"
    
    username, password = SpaceTrack_authenticate()
    siteCred = {'identity': username, 'password': password}

    with requests.Session() as session:
        resp = session.post(uriBase + requestLogin, data=siteCred)
        if resp.status_code != 200:
            raise Exception("POST fail on login")

        satIds = get_satellite_ids(session, uriBase, requestCmdAction, requestFindStarlinks)
        NORAD_ids = [int(i) for i in NORAD_ids]

        for s in NORAD_ids:
            if s in satIds:
                output = get_satellite_data(session, uriBase, requestCmdAction, requestOMMStarlink1, requestOMMStarlink2, s)
                process_satellite_data(output, folder_path, s)
    
    print("Completed session")

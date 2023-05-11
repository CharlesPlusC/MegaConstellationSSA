import os
import configparser
import json
import time
import requests
from datetime import datetime, timedelta

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

def process_satellite_data(output, folder_path, s, raw):
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

            if raw:
                f.write(spacetrack_dict['TLE_LINE1'] + '\n' + spacetrack_dict['TLE_LINE2'] + '\n')
            else:
        # Parse the first line
                    tle_dict["line number"] = line_one[0]
                    tle_dict["satellite catalog number"] = line_one[2:7]
                    tle_dict["classification"] = line_one[7]
                    tle_dict["International Designator(launch year)"] = line_one[9:11]
                    tle_dict["International Designator (launch num)"] = line_one[11:14]
                    tle_dict["International Designator (piece of launch)"] = line_one[
                        14:17
                    ]
                    tle_dict["epoch year"] = line_one[18:20]
                    tle_dict["epoch day"] = line_one[20:32]
                    tle_dict[
                        "first time derivative of mean motion(ballisitc coefficient)"
                    ] = line_one[33:43]
                    tle_dict[
                        "second time derivative of mean motion(delta-dot)"
                    ] = line_one[44:52]
                    tle_dict["bstar drag term"] = line_one[53:61]
                    tle_dict["ephemeris type"] = line_one[62]
                    tle_dict["element number"] = line_one[63:68]
                    tle_dict["checksum"] = line_one[68:69]

                    # Parse the second line (ignore the line number,
                    # satellite catalog number, and checksum)
                    tle_dict["inclination"] = line_two[8:16]
                    tle_dict["right ascension of the ascending node"] = line_two[17:25]
                    tle_dict["eccentricity"] = line_two[26:33]
                    tle_dict["argument of perigee"] = line_two[34:42]
                    tle_dict["mean anomaly"] = line_two[43:51]
                    tle_dict["mean motion"] = line_two[52:63]
                    tle_dict["revolution number at epoch"] = line_two[63:68]

                    #convert epoch to julian day
                    year, decimal_days = tle_dict['epoch year'], tle_dict['epoch day']
                    year = 2000 + int( year)
                    decimal_days = float(decimal_days)
                    #convert the year and decimal days to a datetime object
                    date = datetime(year, 1, 1) + timedelta(days=decimal_days - 1)
                    #convert the datetime object to MJD
                    epoch_mjd = (date - datetime(1858, 11, 17)).total_seconds() / 86400.0
                    #convert to julian date
                    epoch_jd = epoch_mjd + 2400000.5

                    kep_elems = tle_convert(tle_dict)

                    x_car, y_car, z_car, u_car, v_car, w_car = kep2car(
                        a=kep_elems["a"],
                        e=kep_elems["e"],
                        i=kep_elems["i"],
                        w=kep_elems["RAAN"],
                        W=kep_elems["arg_p"],
                        V=kep_elems["true_anomaly"],
                    )
                    state_i = [x_car, y_car, z_car, u_car, v_car, w_car]
                    
                    f.write(
                        "sat:"+ str(tle_dict["satellite catalog number"])+ "ECI:"+ str(state_i)+ "jd_time:"+ str(epoch_jd)+ "\n")
            f.close()# Process the TLE data and write it to the file
    session.close()
    print("Completed Session")

def NORAD_TLE_History(NORAD_ids, constellation, folder_path="/Users/charlesc/Documents/GitHub/Astrodynamics/data/TLEs_Raw/", raw=False):
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
                process_satellite_data(output, folder_path, s, raw)
    
    print("Completed session")

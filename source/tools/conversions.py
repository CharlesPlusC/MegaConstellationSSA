"""
Time and coordinate conversions
"""
import numpy as np
import warnings
from astropy import units as u
import astropy.units as units
from astropy.time import Time
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import datetime
from astropy.coordinates import GCRS, ITRS, CartesianRepresentation, CartesianDifferential, SkyCoord, GCRS, CIRS, TEME, TETE, ITRS, ICRS
from pyproj import Transformer
from typing import List, Tuple

def kep2car(a, e, i, w, W, V):
    # Suppress the UserWarning for true anomaly wrapping
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        
        # Create an Orbit object from the Keplerian elements
        orbit = Orbit.from_classical(Earth,
                                     a * u.km,
                                     e * u.one,
                                     i * u.rad,
                                     w * u.rad,
                                     W * u.rad,
                                     V * u.rad,
                                     epoch=Time.now())

    # Get the position and velocity vectors in ECI frame
    pos_vec = orbit.r.value
    vel_vec = orbit.v.value

    # Extract the Cartesian coordinates and velocities
    x, y, z = pos_vec
    vx, vy, vz = vel_vec

    return x, y, z, vx, vy, vz

def eci2ecef_astropy(eci_pos, eci_vel, mjd):
    # Convert MJD to isot format for Astropy
    time_utc = Time(mjd, format="mjd", scale='utc')

    # Convert ECI position and velocity to ECEF coordinates using Astropy
    eci_cartesian = CartesianRepresentation(eci_pos.T * u.km)
    eci_velocity = CartesianDifferential(eci_vel.T * u.km / u.s)
    gcrs_coords = GCRS(eci_cartesian.with_differentials(eci_velocity), obstime=time_utc)
    itrs_coords = gcrs_coords.transform_to(ITRS(obstime=time_utc))

    # Get ECEF position and velocity from Astropy coordinates
    ecef_pos = np.column_stack((itrs_coords.x.value, itrs_coords.y.value, itrs_coords.z.value))
    ecef_vel = np.column_stack((itrs_coords.v_x.value, itrs_coords.v_y.value, itrs_coords.v_z.value))

    return ecef_pos, ecef_vel

def ecef_to_lla(x, y, z):
    # Convert input coordinates to meters
    x_m, y_m, z_m = x * 1000, y * 1000, z * 1000
    
    # Create a transformer for converting between ECEF and LLA
    transformer = Transformer.from_crs(
        "EPSG:4978", # WGS-84 (ECEF)
        "EPSG:4326", # WGS-84 (LLA)
        always_xy=True # Specify to always return (X, Y, Z) ordering
    )

    # Convert coordinates
    lon, lat, alt_m = transformer.transform(x_m, y_m, z_m)

    # Convert altitude to kilometers
    alt_km = alt_m / 1000

    return lat, lon, alt_km

def jd_to_utc(jd):
    """Converts Julian Date to UTC time tag(datetime object) using Astropy"""
    #convert jd to astropy time object
    time = Time(jd, format='jd', scale='utc')
    #convert astropy time object to datetime object
    utc = time.datetime
    return utc

def utc_jd_date(day, month, year, hours,minutes,seconds, mjd = False, midnight = False):
    """Given a day, month, year, hours, and seconds, return the Julian Date of that time"""

    if midnight == True:
    #convert to datetime object (wihtout hours and seconds)
        date = datetime.datetime(year, month, day)
    else:
        #convert to datetime object (wiht hours and seconds)
        date = datetime.datetime(year, month, day, hours,minutes, seconds)

    #convert to mjd
    mjd = (date - datetime.datetime(1858, 11, 17)).total_seconds() / 86400.0
    if mjd == True:
        return mjd
    else:
        #convert to jd
        jd = mjd + 2400000.5
        return jd
    
def midnight_jd_date(day, month, year, mjd=False):
    """Given a day, month, and year, return the Julian Date of the midnight of that day
    
    Example: today_midnight_jd = midnight_jd_date(11,11,2022)
    """
    #convert to datetime object
    date = datetime.datetime(year, month, day)

    #convert to mjd
    mjd = (date - datetime.datetime(1858, 11, 17)).total_seconds() / 86400.0
    if mjd == True:
        return mjd
    else:
        #convert to jd
        jd = mjd + 2400000.5
        return jd
    
def HCL_diff(eph1,eph2):
    
    """
    Calculate the Height, Cross-Track and Along-Track differences at each time step between two ephemerides.

    Args:
        eph1 (array): list or array of state vectors for a satellite
        eph2 (array): list or array of state vectors for another satellite

    Returns:
        arrays: Three arrays, each containing the height, cross-track and along-track differences at each time step.
    """
    #check that the starting conditions are the same
    # if (eph1[0][0:3]) != (eph2[0][0:3]) or (eph1[0][3:6]) != (eph2[0][3:6]):
    #     warnings.warn('The two orbits do not have the same starting conditions. Make sure this is intentional.')

    H_diffs = []
    C_diffs = []
    L_diffs = []

    for i in range(0, len(eph1), 1):
        #calculate the HCL difference at each time step
        
        r1 = np.array(eph1[i][0:3])
        r2 = np.array(eph2[i][0:3])
        
        v1 = np.array(eph1[i][3:6])
        v2 = np.array(eph2[i][3:6])
        
        unit_radial = r1/np.linalg.norm(r1)
        unit_cross_track = np.cross(r1, v1)/np.linalg.norm(np.cross(r1, v1))
        unit_along_track = np.cross(unit_radial, unit_cross_track)

        #put the three unit vectors into a matrix
        unit_vectors = np.array([unit_radial, unit_cross_track, unit_along_track])

        #subtract the two position vectors
        r_diff = r1 - r2

        #relative position in HCL frame
        r_diff_HCL = np.matmul(unit_vectors, r_diff)

        #height, cross track and along track differences
        h_diff = r_diff_HCL[0]
        c_diff = r_diff_HCL[1]
        l_diff = r_diff_HCL[2]

        H_diffs.append(h_diff)
        C_diffs.append(c_diff)
        L_diffs.append(l_diff)

    return H_diffs, C_diffs, L_diffs

def alt_series(ephemeris):
    """Given a list of positions, return a list of altitudes. Input must be a n*3 array

    Args:
        ephemeris (list): list of positions

    Returns:
        list: list of altitudes
    """
    # Assuming spherical Earth
    Re = 6378.137 #km
    
    #check input in correct format
    if len(ephemeris[0]) != 3:
        print('Input must be a n*3 array')
        return

    alts = []
    for i in range(0, len(ephemeris), 1):
        r_vec = ephemeris[i][1] #position vector
        r_mag = np.abs(np.linalg.norm(r_vec)) #magnitude of position vector
        alt = r_mag-Re #altitude
        alts.append(alt) #add altitude to list

    return alts

def dist_3d(state1, state2):
    """Calculates the distance between two 3D vectors.

    Args:
        state1 (array-like): first position vector (x,y,z)
        state2 (array-like): second position vector (x,y,z)

    Returns:
        float: cartesian distance between the two vectors.
    """
    x1, y1, z1 = state1[0], state1[1], state1[2]
    x2, y2, z2 = state2[0], state2[1], state2[2]
    dist = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    return dist


# Create a transformer for converting between ECEF and LLA
transformer = Transformer.from_crs(
    "EPSG:4978",  # WGS-84 (ECEF)
    "EPSG:4326",  # WGS-84 (LLA)
    always_xy=True  # Specify to always return (X, Y, Z) ordering
)

def ecef_to_lla(x: List[float], y: List[float], z: List[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert Earth-Centered, Earth-Fixed (ECEF) coordinates to Latitude, Longitude, Altitude (LLA).
    
    Args:
        x: x coordinates in km.
        y: y coordinates in km.
        z: z coordinates in km.
    
    Returns:
        lat: Latitudes in degrees.
        lon: Longitudes in degrees.
        alt_km: Altitudes in km.
    """
    # Convert input coordinates to meters
    x_m = np.array(x) * 1000
    y_m = np.array(y) * 1000
    z_m = np.array(z) * 1000   

    # Convert coordinates
    lon, lat, alt_m = transformer.transform(x_m, y_m, z_m)

    # Convert altitude to kilometers
    alt_km = np.array(alt_m) / 1000

    return lat, lon, alt_km

def eci2latlon(eci_positions: List[List[float]], eci_velocities: List[List[float]], mjd_times: List[float]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert Earth-Centered Inertial (ECI) coordinates to Latitude and Longitude.
    
    Args:
        eci_positions: List of ECI positions.
        eci_velocities: List of ECI velocities.
        mjd_times: List of modified Julian dates.
    
    Returns:
        lats: Latitudes in degrees.
        lons: Longitudes in degrees.
    """
    # Check that the first dimensions of eci_positions and eci_velocities are the same
    if len(eci_positions) != len(eci_velocities):
        raise ValueError('eci_positions and eci_velocities must have the same first dimension')

    # Check that each eci_position and eci_velocity is of length 3
    for eci_pos, eci_vel in zip(eci_positions, eci_velocities):
        if len(eci_pos) != 3 or len(eci_vel) != 3:
            raise ValueError('Each eci_pos and eci_vel must be of length 3')

    ecef_positions, _ = eci2ecef_astropy(np.array(eci_positions), np.array(eci_velocities), np.array(mjd_times))
    
    # Convert all ecef_positions to lla_coords at once
    x = [pos[0] for pos in ecef_positions]
    y = [pos[1] for pos in ecef_positions]
    z = [pos[2] for pos in ecef_positions]

    lla_coords = np.array(ecef_to_lla(x, y, z))

    lats = lla_coords[0]
    lons = lla_coords[1]

    return lats, lons

def doy_to_dom_month(year, doy):
    d = datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1)
    day_of_month = d.day
    month = d.month
    return day_of_month, month

def parse_spacex_datetime_stamps(timestamps):
    """Parse SpaceX ephemeris datetime stamps into year, day of year, hour, minute, second."""
    
    # make an array where we will store the year, day of year, hour, minute, and second for each timestamp
    parsed_tstamps = np.zeros((len(timestamps), 7))
    
    for i in range(0, len(timestamps), 1):
        tstamp_str = str(timestamps[i])
        # year is first 4 digits
        year = tstamp_str[0:4]
        # day of year is next 3 digits
        dayofyear = tstamp_str[4:7]
        #convert day of year to day of month and month number
        day_of_month, month = doy_to_dom_month(int(year), int(dayofyear))
        # hour is next 2 digits
        hour = tstamp_str[7:9]
        # minute is next 2 digits
        minute = tstamp_str[9:11]
        # second is next 2 digits
        second = tstamp_str[11:13]
        # milisecond is next 3 digits
        milisecond = tstamp_str[14:16]
        # add the parsed timestamp to the array
        parsed_tstamps[i] = ([int(year), int(month), int(day_of_month), int(hour), int(minute), int(second), int(milisecond)])

    return parsed_tstamps

def yyyy_mm_dd_hh_mm_ss_to_jd(year, month, day, hour, minute, second, milisecond):
    """Convert year, month, day, hour, minute, second to datetime object and then to julian date."""
    dt_obj = datetime.datetime(year, month, day, hour, minute, second, milisecond*1000)
    jd = Time(dt_obj).jd
    return jd

def TEME_to_MEME(x, y, z, u, v, w, jd_time):
    """ convert from the ECI frame used for the NORAD two-line elements: sometimes called true equator, mean equinox (TEME) although it does not use the conventional mean equinox to the mean equator and mean equinox (MEME i.e. GCRS) frame used by the spacex ephemeris.
    """
    # convert to astropy time object
    astropy_time = Time(jd_time, format='jd')
    # convert to astropy skycoord object
    skycoord = SkyCoord(x, y, z, unit='km', representation_type='cartesian', frame=TEME(obstime=astropy_time))
    # convert to GCRS frame
    gcrs = skycoord.transform_to(GCRS(obstime=astropy_time))
    # convert to cartesian coordinates
    x, y, z = gcrs.cartesian.xyz.to(units.km)
    # convert to astropy skycoord object
    skycoord = SkyCoord(u, v, w, unit='km/s', representation_type='cartesian', frame=TEME(obstime=astropy_time))
    # convert to GCRS frame
    gcrs = skycoord.transform_to(GCRS(obstime=astropy_time))
    # convert to cartesian coordinates
    u, v, w = gcrs.cartesian.xyz.to(units.km/units.s)
    
    return x.value, y.value, z.value, u.value, v.value, w.value
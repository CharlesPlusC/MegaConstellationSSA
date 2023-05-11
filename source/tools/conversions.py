"""
Time and coordinate conversions
"""

import numpy as np
import warnings
from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import datetime

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
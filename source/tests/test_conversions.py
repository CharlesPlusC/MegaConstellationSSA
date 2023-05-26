"""Pytest functions."""
import numpy as np
import pytest
import datetime
from ..tools.conversions import HCL_diff, alt_series, dist_3d, ecef_to_lla, eci2latlon, doy_to_dom_month, \
                                jd_to_mjd, parse_spacex_datetime_stamps, yyyy_mm_dd_hh_mm_ss_to_jd, kep2car, \
                                car2kep, eci2ecef_astropy, jd_to_utc

def test_HCL_diff():
    """Test the HCL_diff function."""
    eph1 = np.array([[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]])
    eph2 = np.array([[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]])
    H_diffs, C_diffs, L_diffs = HCL_diff(eph1, eph2)
    assert H_diffs == [0, 0]
    assert C_diffs == [0, 0]
    assert L_diffs == [0, 0]

def test_alt_series():
    """Test the alt_series function."""
    ephemeris = [[1, [1, 1, 6378.137 + 10]], [2, [1, 1, 6378.137 + 20]]]
    assert alt_series(ephemeris) == [10, 20]

def test_dist_3d():
    """Test the dist_3d function."""
    state1 = [1, 2, 3]
    state2 = [4, 5, 6]
    assert dist_3d(state1, state2) == np.sqrt((1-4)**2 + (2-5)**2 + (3-6)**2)

def test_ecef_to_lla():
    """Test the ecef_to_lla function."""
    x = [1, 2, 3]
    y = [4, 5, 6]
    z = [7, 8, 9]
    lat, lon, alt_km = ecef_to_lla(x, y, z)
    assert isinstance(lat, np.ndarray)
    assert isinstance(lon, np.ndarray)
    assert isinstance(alt_km, np.ndarray)

def test_eci2latlon():
    """Test the eci2latlon function."""
    eci_positions = [[1, 2, 3], [4, 5, 6]]
    eci_velocities = [[7, 8, 9], [10, 11, 12]]
    mjd_times = [13, 14]
    lats, lons = eci2latlon(eci_positions, eci_velocities, mjd_times)
    assert isinstance(lats, np.ndarray)
    assert isinstance(lons, np.ndarray)

def test_doy_to_dom_month():
    """Test the doy_to_dom_month function."""
    assert doy_to_dom_month(2023, 1) == (1, 1)
    assert doy_to_dom_month(2023, 365) == (31, 12)
    assert doy_to_dom_month(2024, 366) == (31, 12)  # Leap year

def test_jd_to_mjd():
    """Test the jd_to_mjd function."""
    # 2451545 JD is 51544.5 MJD
    assert jd_to_mjd([2451545]) == [51544.5]

def test_parse_spacex_datetime_stamps():
    """Test the parse_spacex_datetime_stamps function."""
    # Assuming timestamp in the format YYYYDDDHHMMSS
    timestamp = ['2022014123456']
    expected = np.array([[2022, 1, 14, 12, 34, 56, 0]])
    np.testing.assert_array_equal(parse_spacex_datetime_stamps(timestamp), expected)

def test_yyyy_mm_dd_hh_mm_ss_to_jd():
    """Test the yyyy_mm_dd_hh_mm_ss_to_jd function."""
    # Check Julian date for a known date
    jd = yyyy_mm_dd_hh_mm_ss_to_jd(2000, 1, 1, 12, 0, 0, 0)
    # For this date and time, JD is 2451545.0
    assert jd == 2451545.0

def test_kep2car():
    # Test using a known orbit
    a = 42164.0  # km, geostationary orbit
    e = 0.0
    i = 0.0  # radians
    w = 0.0  # radians
    W = 0.0  # radians
    V = 0.0  # radians

    x, y, z, vx, vy, vz = kep2car(a, e, i, w, W, V)
    assert pytest.approx(x, abs=1e-9) == -42164.0
    assert pytest.approx(y, abs=1e-9) == 0.0
    assert pytest.approx(z, abs=1e-9) == 0.0

def test_car2kep():
    # Test using a known orbit
    x = -42164.0  # km, geostationary orbit
    y = 0.0
    z = 0.0
    vx = 0.0
    vy = -3.074  # km/s, geostationary orbit
    vz = 0.0

    a, e, i, w, W, V = car2kep(x, y, z, vx, vy, vz, deg=False)
    assert pytest.approx(a, abs=1e-9) == 42164.0
    assert pytest.approx(e, abs=1e-9) == 0.0
    assert pytest.approx(i, abs=1e-9) == 0.0

def test_eci2ecef_astropy():
    # Test using a known position and velocity
    eci_pos = np.array([[42164.0, 0.0, 0.0]])  # km
    eci_vel = np.array([[0.0, 3.074, 0.0]])  # km/s
    mjd = 51544.5  # 2000-01-01 00:00:00 UTC

    ecef_pos, ecef_vel = eci2ecef_astropy(eci_pos, eci_vel, mjd)
    # Just testing x coordinate here, should test all coordinates and velocity
    assert pytest.approx(ecef_pos[0, 0], abs=1e-9) == 42164.0

def test_ecef_to_lla():
    # Test using a known position
    x = [42164.0]
    y = [0.0]
    z = [0.0]

    lat, lon, alt = ecef_to_lla(x, y, z)
    assert pytest.approx(lat[0], abs=1e-9) == 0.0
    assert pytest.approx(lon[0], abs=1e-9) == 0.0
    assert pytest.approx(alt[0], abs=1e-9) == 35786.0  # Geostationary altitude

def test_jd_to_utc():
    # Test using a known date
    jd = 2451544.5  # 2000-01-01 00:00:00 UTC
    utc = jd_to_utc(jd)
    assert utc == datetime.datetime(2000, 1, 1, 0, 0, 0)
"""Script to download relevant TLE data from Space-Track."""
from tools.tletools import download_tle_history, load_satellite_lists, NORAD_list_update

def download_data_for_analysis():
    """Download the TLEs from Space-track for the NORAD IDs in the satellite_lists.json file."""
    # load the satellite lists from the source/satellite_lists.json file
    satellite_lists = load_satellite_lists()

    # Access the lists using the satellite_lists dictionary
    norad_OW_L5 = satellite_lists["oneweb"]["L5"]
    norad_OW_L4 = satellite_lists["oneweb"]["L4"]
    norad_OW_L6 = satellite_lists["oneweb"]["L6"]

    norad_SL_L28 = satellite_lists["starlink"]["L28"]
    norad_SL_L30 = satellite_lists["starlink"]["L30"]
    norad_SL_L36 = satellite_lists["starlink"]["L36"]

    #merge the lists above into constellation-wide lists
    SL_norads = norad_SL_L28 + norad_SL_L36 + norad_SL_L30
    OW_norads = norad_OW_L4 + norad_OW_L5 + norad_OW_L6

    # download the list of NORAD IDs for the specified constellations
    NORAD_list_update('oneweb')
    NORAD_list_update('starlink')

    # fetch the available TLE archive for the specified NORAD IDs
    download_tle_history(SL_norads, 'starlink') 
    download_tle_history(OW_norads, 'oneweb')

if __name__ == "__main__":
    download_data_for_analysis()
import os
import sys
from tools.analysis_tools import NORAD_vs_SUP_TLE_analysis, TLE_analysis_to_df
from tools.tletools import load_satellite_lists
import cProfile

def main ():
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
    
    #convert the list of ints to strings for the function
    SL_norads = [str(i) for i in SL_norads]
    OW_norads = [str(i) for i in OW_norads]

    all_norads = SL_norads + OW_norads
    print("Norad IDS specified: ", all_norads)
    # Unhash to run the analysis
    # NORAD_vs_SUP_TLE_analysis(NORADS = all_norads)

    Oneweb_dfs, Starlink_dfs = TLE_analysis_to_df()

    print("columns: ", Oneweb_dfs[0].columns)

if __name__ == "__main__":
    main()
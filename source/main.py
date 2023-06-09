"""Main Analysis Script."""
import os
import sys
from tools.analysis_tools import TLE_rate_dicts, NORAD_vs_SUP_TLE_analysis, TLE_analysis_to_df, launch_specific_stats, sup_gp_op_benchmark, TLE_arglat_dict
from tools.tletools import load_satellite_lists
from tools.plotting_tools import plot_map_diffs_smallvals_subplot, plot_tle_rate_analysis, plot_arglat_analysis, plot_altitude_timeseries, plot_fft_comparison, plot_diff_subplots, plot_diff_hist, plot_launch_latlon_diffs, plot_ground_tracks, plot_map_diffs_smallvals_all, benchmark_plot
import cProfile
import numpy as np

def main ():
    """Run all the analyses available."""
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
    # Unhash to run the analysis
    # NORAD_vs_SUP_TLE_analysis(NORADS = all_norads) #TODO:specifying the NOARDS doesnt actually seem to do anytbing...

    #load in the analysis data
    Oneweb_dfs, Starlink_dfs = TLE_analysis_to_df(NORAD_IDs = all_norads) # just taking the first two dataframes for speed

    # # plot the altitude time series for the selected list of dataframes
    # plot_altitude_timeseries(Starlink_dfs, show=True)
    # plot_altitude_timeseries(Oneweb_dfs, show=True)

    all_dfs = Oneweb_dfs + Starlink_dfs # combine the two lists of dataframes

    # # calculate the launch-specific stats and store them in .csv file
    # launch_specific_stats(all_dfs)

    # # plot the FFTs for the selected list of dataframes
    # plot_fft_comparison(all_dfs, show=True)

    # # Histogram of the difference between the NORAD and SUP TLEs
    # plot_diff_hist(all_dfs, show=False)

    # # Error as a function of geographic location (lat/lon)
    # plot_launch_latlon_diffs(all_dfs, show=False, criteria=1) #errors within 1 SD from the mean 
    # plot_launch_latlon_diffs(all_dfs, show=False, criteria=2) #errors within 2 SD from the mean 

    # # Plot the ground tracks for the selected list of dataframes
    # plot_ground_tracks(all_dfs, show=True)

    # # Plot the difference between the NORAD and SUP TLEs projected onto the ground track
    # plot_map_diffs_smallvals_all(all_dfs, show=False, criteria=1) #errors within 1 SD from the mean

    plot_map_diffs_smallvals_subplot(all_dfs, show=False, criteria=1) #errors within 1 SD from the mean
    
    # # Histogram of the argument of latitude at which each TLE is produced for all available NORAD IDs as a function of source
    # plot_arglat_analysis()

    # # Histogram of the TLE latency for each constellation
    # plot_tle_rate_analysis()

    # # Benchmarking the NORAD, Suppelemental, and Operator ephemerides for 3 Starlink satellites. Uses the data in external/ephem_TLE_compare
    # benchmark_plot() 

    # TODO if time: TLE age vs error
    # TODO if time: H/C/L correlation
    # TODO if time: Solar Flux vs error

if __name__ == "__main__":
    main()
    
#   Unhash to run the profiler
    # cProfile.run('main()')
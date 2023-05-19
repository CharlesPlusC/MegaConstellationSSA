"""Set of plotting tools for the project. Most of these are based on the use of lists of pandas dataframes."""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
import json
from scipy import signal
from typing import Dict, List

#local imports
from .analysis_tools import compute_fft

# Set the default font size for the plots
mpl.rcParams['font.size'] = 11

# Dictionary of colours for the different launches to be consistent across plots
launch_colour_dict = {'L4': 'xkcd:blue', 'L28': 'xkcd:dark red', 'L5': 'xkcd:azure', 'L36': 'xkcd:orange', 'L6': 'xkcd:light blue', 'L30': 'xkcd:coral'}
# list of the different types of differences to be plotted
diff_types = ['h_diffs', 'c_diffs', 'l_diffs', 'cart_pos_diffs'] 
              
def plot_altitude_timeseries(dfs, json_filepath='external/selected_satellites.json', show=False):
    """
    Plots altitude time series for different satellites from given dataframes.

    Parameters:
    dfs (list): A list of pandas DataFrames, each containing data for a different satellite.
    json_filepath (str): The path to the JSON file containing the selected satellite numbers.
    """

    with open(json_filepath, 'r') as f:
        selected_satellites = json.load(f)

    sat_nums = []
    for df in dfs:
        sat_nums.append(df['NORAD_ID'].unique()[0])

    eph_alts_sup = []
    times = []
    launch_nos = []
    colours = []
    added_labels = set() # keep track of labels already added to the legend
    constellations = {} # keep track of constellations for each satellite number

    for sat_num in sat_nums:
        for df in dfs:
            if sat_num == df['NORAD_ID'].unique()[0]:
                eph_alts_sup.append(df['eph_alts_sup'])
                times.append(df['times'] - 2400000.5) # convert from JD to MJD for readability

                # Find the constellation and launch number by comparing with the selected satellites
                for constellation, launches in selected_satellites.items():
                    for launch, norad_ids in launches.items():
                        if int(df['NORAD_ID'][0]) in norad_ids:
                            launch_nos.append(launch)
                            colours.append(launch_colour_dict[launch]) # use the launch colour dictionary to get the colour for the launch
                            constellations[sat_num] = constellation

    fig, ax = plt.subplots(figsize=(8, 5))
    for i in range(len(eph_alts_sup)):
        # only add label if it hasn't been added before
        label = launch_nos[i] if launch_nos[i] not in added_labels else ""
        ax.scatter(x=times[i], y=eph_alts_sup[i], color=colours[i], label=label, alpha=0.2, s=0.1)
        added_labels.add(launch_nos[i]) # mark this label as added
    ax.set_title(constellations[sat_nums[0]]) # set title to the constellation of the first satellite number
    ax.set_ylabel('Altitude (Km)')
    ax.set_xlabel('Time (MJD)')
    ax.set_xlim(59390, 59940) # window of interest for this project
    legend = ax.legend(loc='upper right')
    # set the alpha value of the legend markers higher
    for lh in legend.legendHandles: 
        lh.set_alpha(1)
        # and make the marker size bigger
        lh._sizes = [4]
    ax.grid()
    plt.tight_layout()
    plt.savefig(f'output/plots/altitude/altitude_tseries_{constellations[sat_nums[0]]}.png', dpi=400)
    if show:
        plt.show()


def plot_fft_comparison(list_of_dfs: List[pd.DataFrame], diff_types: List[str] = diff_types, launch_colour_dict: Dict[str, str] = launch_colour_dict, show: bool = True) -> None:
    """
    Plots the power spectral density of each dimension (H/C/L/3D) in the differences 
    between NORAD and operator TLEs for each launch.
    
    Args:
        list_of_dfs (List[pd.DataFrame]): List of dataframes each representing a specific launch.
        diff_types (List[str]): List of difference types for which the plot needs to be drawn.
        launch_colour_dict (Dict[str, str]): Dictionary mapping launch id to colour for plotting.
        show (bool, optional): Flag indicating whether to display the plot. Defaults to True.

    Returns:
        None: The function performs plotting operation and does not return any value.
    """
    grouped_dfs = {}
    for df in list_of_dfs:
        constellation = str(df['constellation'][0])
        if constellation not in grouped_dfs:
            grouped_dfs[constellation] = [df]
        else:
            grouped_dfs[constellation].append(df)

    # Now loop through each constellation and plot its dataframes
    for constellation, dfs in grouped_dfs.items():
        for diff_type in diff_types:
            figure, axis = plt.subplots(2, 1)
            figure.set_size_inches(7, 10)

            for df in dfs:
                launch = 'L' + str(df['launch_no'][0])
                col = launch_colour_dict.get(launch, 'black')

                # In the time domain
                axis[0].scatter(df['times'], df[diff_type], alpha=0.1, s=1, color=col)

                # In the frequency domain
                diff_fft, diff_psd = compute_fft(df, diff_type)
                axis[1].plot(diff_fft, diff_psd, alpha=0.2, color=col)

            title_dict = {
                'l_diffs': r'$\Delta$ L: NORAD- and Operator-TLE Derived positions for {}'.format(constellation),
                'cart_pos_diffs': r'$\Delta$ 3D: NORAD- and Operator-TLE Derived positions for {}'.format(constellation),
                'c_diffs': r'$\Delta$ C: NORAD- and Operator-TLE Derived positions for {}'.format(constellation),
                'h_diffs': r'$\Delta$ H: NORAD- and Operator-TLE Derived positions for {}'.format(constellation)
            }

            axis[0].set_ylabel('Delta(Km)')
            axis[1].set_ylabel('Power Spectral Density (dB)')
            axis[1].set_xlabel('Frequency (days^-1)')
            axis[0].set_xlabel('Julian Day')
            axis[1].set_title('Power spectral density of the differences between the TLEs')
            axis[1].grid(True)
            axis[0].grid(True)
            axis[0].set_xlim(np.min(df['times']), np.max(df['times']))
            axis[0].set_title(title_dict.get(diff_type, ''))

            axis[1].set_xticks(np.arange(1, 40, 2))
            axis[1].set_xticklabels(np.arange(1, 40, 2))
            axis[1].set_xlim(0, 40)

            plt.savefig('output/plots/Fourier_analysis/{}_{}_fft.png'.format(constellation, diff_type), dpi=300)
            
            if show == True:
                plt.show()
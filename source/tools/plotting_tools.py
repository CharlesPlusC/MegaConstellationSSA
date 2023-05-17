"""Set of plotting tools for the project. Most of these are based on the use of lists of pandas dataframes."""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
import json
from scipy import signal

mpl.rcParams['font.size'] = 11

# Dictionary of colours for the different launches to be consistent across plots
launch_colour_dict = {'L4': 'xkcd:blue', 'L28': 'xkcd:dark red', 'L5': 'xkcd:azure', 'L36': 'xkcd:orange', 'L6': 'xkcd:light blue', 'L30': 'xkcd:coral'}
# list of the different types of differences to be plotted
diff_types = ['l_diffs', 'c_diffs', 'l_diffs', 'cart_pos_diffs'] 
              
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


def plot_fft_compare(diff_type, launch_data_dict, launch_colour_dict, show = False):
    figure, axis = plt.subplots(2, 1, figsize = (7, 10))
    for launch, data in launch_data_dict.items():
        colour = launch_colour_dict[launch]
        freqs, psd = signal.welch(data[diff_type], fs=1, nperseg=4096)
        axis[0].semilogy(freqs, psd, color=colour)
        axis[1].semilogy(freqs, psd, color=colour)
        axis[1].legend([Patch(facecolor = colour, edgecolor = colour, label = launch) for launch, colour in launch_colour_dict.items()])

    axis[0].set_xticks(np.arange(0, 10, 1))
    axis[0].set_xlim(0, 10)
    axis[0].set_ylim(0, 80)
    axis[0].grid(True)
    
    axis[1].set_xticks(np.arange(1, 40, 2))
    axis[1].set_xlim(0, 40)
    axis[1].set_ylim(0, 80)
    axis[1].grid(True)
    axis[1].set_xlabel('Frequency (1/Day)', fontsize=12)
    
    if diff_type == 'l_diffs':
        figure.suptitle(r'$\Delta$ L: Power Spectral Density Comparison')
    elif diff_type == 'cart_pos_diffs':
        figure.suptitle(r'$\Delta$ 3D: Power Spectral Density Comparison')
    elif diff_type == 'c_diffs':
        figure.suptitle(r'$\Delta$ C: Power Spectral Density Comparison')
    elif diff_type == 'h_diffs':
        figure.suptitle(r'$\Delta$ H: Power Spectral Density Comparison')

    figure.text(0, 0.5, 'Power Spectral Density (dB)', va='center', rotation='vertical')
    figure.subplots_adjust(hspace=0.5)
    plt.tight_layout()

    plt.savefig('output/plots/Fourier_analysis/fft_compare_' + diff_type + '.png', dpi=300)
    if show == True:
        plt.show()

def run_fouier_analysis(diff_types, constellations, launch_data_dict, launch_colour_dict):
    for diff_type in diff_types:
        for constellation in constellations:
            # Assuming that launch_data_dict and launch_colour_dict are dictionaries where the keys are constellations
            constellation_data = launch_data_dict[constellation]
            constellation_colour = launch_colour_dict[constellation]
            
            plot_fft_compare(diff_type, constellation_data, constellation_colour)

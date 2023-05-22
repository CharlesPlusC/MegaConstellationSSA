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
constellation_colour_dict = {'oneweb': 'xkcd:azure', 'starlink': 'xkcd:coral'}
# list of the different types of differences to be plotted
diff_types = ['h_diffs', 'c_diffs', 'l_diffs', 'cart_pos_diffs'] 

# Mapping of satellite numbers to constellations and launch numbers
json_filepath='external/selected_satellites.json'
with open(json_filepath, 'r') as f:
    selected_satellites = json.load(f)
              
def plot_altitude_timeseries(dfs, show=False):
    """
    Plots altitude time series for different satellites from given dataframes.

    Parameters:
    dfs (list): A list of pandas DataFrames, each containing data for a different satellite.
    json_filepath (str): The path to the JSON file containing the selected satellite numbers.
    """

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
    Plots the frequency and power spectral density  of each dimension (H/C/L/3D) time-series in the differences 
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
            
            if show:
                plt.show()

def plot_diff_subplots(sats_dataframe, diffs='all', show=False):
    # Divide dataframes by constellation
    constellation_dict = {}
    for df in sats_dataframe:
        constellation = df['constellation'][0]
        if constellation not in constellation_dict:
            constellation_dict[constellation] = []
        constellation_dict[constellation].append(df)

    if diffs == 'all':
        diffs_types = np.array(diff_types)
    else:
        diffs_types = np.array(diffs)

    for constellation, sats_dataframe in constellation_dict.items():
        for j in range(0, len(diffs_types), 1):
            diff_type = diffs_types[j]
            sats_dataframe.sort(key=lambda x: x['NORAD_ID'].iloc[0])

            fig, axs = plt.subplots(6, 5, figsize=(15, 9))
            fig.suptitle(r'$\Delta$' + diff_type[0].upper()+ 'between NORAD- and Operator-TLE derived SGP4 orbits for ' + constellation + ' Satellites')

            earliest_time = 59400
            latest_time = 59925
            i = 0
            used_colors = set()  # Store the used colors for legend
            for row in axs:
                for ax in row:
                    if i < len(sats_dataframe):
                        df = sats_dataframe[i]
                        launch = df['launch_no'][0]
                        launch_key = 'L' + str(launch)
                        if launch_key in launch_colour_dict:
                            col = launch_colour_dict[launch_key]
                            used_colors.add(launch_key)
                        else:
                            raise Exception('Launch key not found in launch_colour_dict')

                        ax.scatter(df['MJD'], df[diff_type], s=0.1, alpha=0.2, c=col)
                        ax.set_facecolor('xkcd:grey')
                        ax.set_title('NORAD:' + str(df['NORAD_ID'][0]), fontsize=10)

                        min_max_dict = {'h_diffs': (-2, 2), 'c_diffs': (-1, 1), 'l_diffs': (-20, 20),
                                        'cart_pos_diffs': (0, 20)}
                        min, max = min_max_dict.get(diff_type, (0, 0))
                        ax.set_ylim(min, max)
                        fig.text(0.08, 0.5, r'$\Delta$ ' + diff_type[0].upper() + ' (Km)', ha='center', va='center',
                                 rotation='vertical', fontsize=10)

                        if i > 23:
                            ax.set_xticks(np.linspace(earliest_time, latest_time, 4))
                        else:
                            ax.set_xticks([])
                        if i % 5 == 0:
                            ax.set_yticks([min, 0, max])
                        else:
                            ax.set_yticks([])

                        ax.text(0.025, 0.96,
                                'µ:' + str(round((df[diff_type].mean()) * 1000, 3)) + ' m' + '; ' + 'σ: ' + str(
                                    round(df[diff_type].std() * 1000, 3)) + ' m', transform=ax.transAxes, fontsize=10,
                                verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))

                        fig.text(0.5, 0.02, 'Modified Julian Day', ha='center', va='center', fontsize=12)

                        patch_list = [Patch(facecolor=launch_colour_dict[key], edgecolor=launch_colour_dict[key], label=key)
                                      for key in used_colors]  # Change here, only include used colors
                        fig.legend(patch_list, [key for key in used_colors], loc='center right')  # Change here, only include used colors

                        plt.setp(ax.get_xticklabels(), rotation=50, ha="right", rotation_mode="default")
                        i += 1

            plt.subplots_adjust(hspace=0.3, wspace=0.08)
            plt.savefig('output/plots/timseries_subplots/'+diff_type+'_subplots_'+constellation+'.png', dpi=300)
            if show:
                plt.show()

def plot_diff_hist(sats_dataframe_list, diffs='all', show=False):
    fig, axs = plt.subplots(2, 2, figsize=(8, 5))

    # Define info for each subplot
    subplot_info = {
        'h_diffs': (axs[0, 0], r'$\Delta$ H (Km)', 'Count', [-1, 1]),
        'c_diffs': (axs[0, 1], r'$\Delta$ C (Km)', 'Count', [-1, 1]),
        'l_diffs': (axs[1, 0], r'$\Delta$ L (Km)', 'Count', [-20, 20]),
        'cart_pos_diffs': (axs[1, 1], r'$\Delta$ 3D (Km)', 'Count', [0, 25])
    }

    if diffs != 'all':
        diff_types = diffs
    else:
        diff_types = subplot_info.keys()

    for diff in diff_types:
        ax, xlabel, ylabel, xlims = subplot_info[diff]
        data_dict = {df['constellation'].iloc[0]: df[diff] for df in sats_dataframe_list}

        text_y_offset = 0.15
        text_y = 0.95

        for constellation, data in data_dict.items():
            # Generate histogram with automatic bins
            n, bins, patches = ax.hist(data, bins='auto', alpha=0.5, label=constellation, range=xlims, color=constellation_colour_dict[constellation.lower()])

            # Calculate statistics
            mean, std = np.mean(data), np.std(data)

            # Add statistics to the plot as text
            stats_text = "\n".join([
                f'Mean {constellation}: {mean:.2f}km',
                f'Std {constellation}: {std:.2f}km'
            ])
            ax.text(0.05, text_y, stats_text, transform=ax.transAxes, verticalalignment='top', fontsize=10)
            text_y -= text_y_offset

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, which='major')
        ax.set_yscale('log')
        ax.set_ylim(bottom=1)  # set y lower limit to 1
        ax.set_xlim(xlims)  # set x limits
        ax.legend()

    plt.tight_layout()
    plt.savefig('output/plots/histograms/'+diff+'_hist.png', dpi=300)
    if show:
        plt.show()

def plot_launch_latlon_diffs(sats_dataframe_list: List[pd.DataFrame] = [], show=False, criteria=1):
    
    latslons = ['lats', 'lons']

    # Determine the unique constellations and launches from the JSON data
    constellations = list(selected_satellites.keys())
    launches = list(set([launch for const in constellations for launch in selected_satellites[const].keys()]))

    # Create a subplots grid
    fig, axs = plt.subplots(len(constellations), len(launches), figsize=(5*len(launches), 5*len(constellations)))
    
    # One row per constellation
    no_rows = len(constellations)
    # number of cols is the maximum number of launches in a constellation
    no_cols = max([len(selected_satellites[const].keys()) for const in constellations])
    
    # keep the extra plots blank if there are less launches for one of the constellation

    
    for latlon in latslons:
        for diff in diff_types:
            for df in sats_dataframe_list:
                # Apply criteria
                mean = np.mean(df[diff])
                std = np.std(df[diff])
                df_filtered = df[(df[diff] > mean - criteria*std) & (df[diff] < mean + criteria*std)]
                
                # Assign a color based on the launch
                launch = 'L' + str(df['launch_no'][0])
                constellation = df['constellation'][0]
                color = launch_colour_dict[launch]
                
                # Only plot the data if this constellation-launch combination exists in the JSON data
                if launch in selected_satellites.get(constellation, {}):
                    # Get the correct axis for this constellation and launch
                    row = constellation_to_axis_row[constellation]
                    col = launch_to_axis_col[launch]
                    axis = axs[row, col]

                    # Plot the data
                    axis.scatter(df_filtered[latlon], df_filtered[diff], alpha=0.3, s=0.8, color=color)
                    axis.set_title(f"Δ {diff} vs {latlon.capitalize()}, {criteria}σ")
                    axis.set_xlabel(f'{latlon.capitalize()} (Degrees)')
                    axis.grid(True)

    plt.tight_layout()
    plt.savefig('output/plots/latlon/'+diff+ '_' + latlon + '_diffs.png', dpi=300)
    
    if show:
        plt.show()

if __name__ == "__main__":
    pass
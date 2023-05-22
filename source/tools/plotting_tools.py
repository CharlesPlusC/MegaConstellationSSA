"""Set of plotting tools for the project. Most of these are based on the use of lists of pandas dataframes."""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import json
from scipy import signal
from typing import Dict, List, Union

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
              
def plot_altitude_timeseries(dfs: List[pd.DataFrame], show: bool = False) -> None:
    """
    Generate a time series plot of satellite altitudes.

    This function generates a scatter plot of satellite altitudes as a function of time, with each 
    satellite represented as a different series on the plot. The plot is saved to a specific directory
    and optionally displayed. 

    Parameters
    ----------
    dfs : list of pandas.DataFrame
        A list of dataframes, each containing satellite data for a unique satellite. 
        The columns should include 'NORAD_ID', 'eph_alts_sup', and 'times'.
    show : bool, optional
        If True, the plot is displayed; otherwise, it is only saved. Defaults to False.

    Returns
    -------
    None.

    Notes
    -----
    The function saves the generated plot to the 'output/plots/altitude/' directory and additionally 
    displays it if `show` is set to True.
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
    Generate a plot of the power spectral density (Fourier analysis) of the differences between NORAD 
    and operator TLEs for each launch.

    This function generates a plot in both time and frequency domains for the specified difference types. 
    Each launch is represented by a different series on the plot. The plots are saved to a specific directory
    and optionally displayed. 

    Parameters
    ----------
    list_of_dfs : list of pandas.DataFrame
        A list of dataframes, each containing data for a specific launch. 
        The columns should include 'times', 'constellation', and the difference types.
    diff_types : list of str, optional
        The types of differences to be plotted. Defaults to the global `diff_types`.
    launch_colour_dict : dict, optional
        A dictionary mapping launch ids to colors for plotting. Defaults to the global `launch_colour_dict`.
    show : bool, optional
        If True, the plots are displayed; otherwise, they are only saved. Defaults to True.

    Returns
    -------
    None.

    Notes
    -----
    The function saves the generated plots to the 'output/plots/Fourier_analysis/' directory and additionally 
    displays them if `show` is set to True.
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

def plot_diff_subplots(sats_dataframe: List[pd.DataFrame],
                       diffs: Union[str, List[str]] = 'all',
                       show: bool = True) -> None:
    """
    Generate subplots for specified differences in satellite data per constellation.

    This function generates a grid of subplots for specified differences, or all available differences if 'all' 
    is specified. Each satellite's difference data is plotted in a separate subplot, with the satellites grouped 
    by constellation. 

    Parameters
    ----------
    sats_dataframe : list of pandas.DataFrame
        A list of dataframes, each containing satellite data. 
    diffs : str or list of str, optional
        The differences to plot subplots for. If 'all', plots subplots for all available differences. 
        Otherwise, it should be a list of keys corresponding to the difference columns in the dataframes. 
        Defaults to 'all'.
    show : bool, optional
        If True, the plot is displayed; otherwise, it is only saved. Defaults to True.

    Returns
    -------
    None.

    Raises
    ------
    Exception:
        If a launch key from the dataframes is not found in the 'launch_colour_dict'.

    Notes
    -----
    The function saves the generated plot(s) to the 'output/plots/timeseries_subplots/' directory and additionally 
    displays them if `show` is set to True. Each subplot includes statistical data (mean and standard deviation) 
    presented as text and uses color to distinguish different satellite launches.
    """

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

def plot_diff_hist(sats_dataframe_list: List[pd.DataFrame],
                   diffs: Union[str, List[str]] = 'all',
                   show: bool = False) -> None:
    """
    Generate histograms for specified differences in satellite data.

    This function generates histograms for a list of specified differences or all available differences if 'all'
    is specified. Each difference is plotted in a separate subplot, and the data is segmented by the satellite 
    constellation. The histogram also includes statistical data (mean and standard deviation) presented as text
    in each subplot.

    Parameters
    ----------
    sats_dataframe_list : list of pandas.DataFrame
        A list of dataframes, each containing satellite data, including differences in height (h_diffs), 
        cross-track (c_diffs), along-track (l_diffs), and 3D Cartesian position (cart_pos_diffs).
    diffs : str or list of str, optional
        The differences to plot histograms for. If 'all', plots histograms for all available differences. 
        Otherwise, it should be a list of keys corresponding to the difference columns in the dataframes. 
        Defaults to 'all'.
    show : bool, optional
        If True, the plot is displayed, otherwise it is only saved. Defaults to False.

    Returns
    -------
    None.

    Notes
    -----
    The function saves the generated plot(s) to the 'output/plots/histograms/' directory and additionally displays
    them if `show` is set to True.

    Raises
    ------
    KeyError:
        If the keys of `constellation_colour_dict` do not align with the constellation identifiers in the dataframes.
    """
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
    """
    Generate scatter plots of latitude and longitude differences for satellites from various launches. 

    This function generates a subplot for each unique launch found in the dataframes within `sats_dataframe_list`, 
    plotting the differences in latitude and longitude for each satellite in that launch. The differences plotted 
    are filtered based on a certain criterion, only including values within `criteria` standard deviations of the mean.

    Parameters
    ----------
    sats_dataframe_list : list of pandas.DataFrame, optional
        List of dataframes, each containing satellite information including 'lats', 'lons', 'launch_no', 
        and 'constellation'. Defaults to an empty list.
    show : bool, optional
        If True, the plot is displayed, otherwise it is only saved. Defaults to False.
    criteria : int, optional
        Number of standard deviations from the mean within which to include difference values. Defaults to 1.

    Returns
    -------
    None.

    Notes
    -----
    The function saves the generated plot(s) to the 'output/plots/latlon/' directory and additionally displays
    them if `show` is set to True. The subplot layout is currently specific to a 2x3 configuration.

    Raises
    ------
    KeyError:
        If the keys of `launch_colour_dict` do not align with the launch identifiers in the dataframes.

    """

    #NOTE: this is specific to a 2*3 subplot layout
    #TODO: make this more general

    latslons = ['lats', 'lons']

    # Create a figure with as many subplots as there are unique launches
    launches = ['L' + str(df['launch_no'][0]) for df in sats_dataframe_list]
    unique_launches = sorted(set(launches)) # sorted to ensure consistent mapping
    constellations = sorted(set(df['constellation'][0] for df in sats_dataframe_list)) # getting unique constellations
    fig, axs = plt.subplots(2, len(unique_launches) // 2, figsize=(5 * len(unique_launches) // 2, 5 * 2))

    # Create a dictionary mapping each unique launch to an axis
    launch_to_axis = {launch: axs[i // 3, i % 3] for i, launch in enumerate(unique_launches)}

    for latlon in latslons:
        for diff in diff_types:
            for df in sats_dataframe_list:
                # Apply criteria
                mean = np.mean(df[diff])
                std = np.std(df[diff])
                df_filtered = df[(df[diff] > mean - criteria*std) & (df[diff] < mean + criteria*std)] # remove error more than 1 SD from the mean

                # Assign a color based on the launch
                launch = 'L' + str(df['launch_no'][0])
                color = launch_colour_dict[launch]

                # Create subplot title
                title = f"{df['constellation'][0]} {launch}"

                # Get the correct axis for this launch
                axis = launch_to_axis[launch]

                # Plot the data
                axis.scatter(df_filtered[latlon], df_filtered[diff], alpha=0.3, s=0.8, color=color)
                axis.set_title(title)
                axis.set_xlabel(f'{latlon.capitalize()} (Degrees)')
                axis.grid(True)
                
                # create title for the whole figure
                fig.suptitle(f'{diff.capitalize()} for {latlon.capitalize()}\n {criteria} SD from the mean')

                plt.tight_layout()
            plt.savefig('output/plots/latlon/' + str(diff) + '_' + str(latlon) + '_diffs' + str(criteria) + '_sd' + '.png', dpi=300)

    if show:
        plt.show()

def plot_ground_tracks(list_of_dfs: List[pd.DataFrame] = [], show: bool = False):
    """
    Generate a plot of satellite ground tracks for selected launches from different constellations. 

    The function reads in a list of dataframes, each representing a satellite with its geographic coordinates 
    over time. It then plots the ground track for each satellite on a Basemap plot. The dataframes in 
    `list_of_dfs` should contain columns 'lons', 'lats' and 'launch_no'. The color of each satellite's ground 
    track is determined by its associated launch, as defined in the `launch_colour_dict`. The launch information 
    is loaded from an external JSON file.

    Parameters
    ----------
    list_of_dfs : list of pandas.DataFrame, optional
        List of dataframes each containing satellite information including 'lons', 'lats' and 'launch_no'. 
        Defaults to an empty list.
    show : bool, optional
        If True, the plot is displayed, otherwise it is only saved. Defaults to False.

    Returns
    -------
    None.

    Notes
    -----
    The function saves the generated plot to 'output/plots/ground_tracks/g_tracks.png' and additionally 
    displays it if `show` is set to True.
    
    Raises
    ------
    KeyError:
        If the keys of `launch_colour_dict` do not align with the launch identifiers in the JSON file.

    """
    # Load selected satellites from JSON file
    json_filepath = 'external/selected_satellites.json'
    with open(json_filepath, 'r') as f:
        selected_satellites = json.load(f)
    
    # Prepare the figure
    fig = plt.figure(figsize=(8,5))
    plt.rcParams['font.size'] = 11.0
    m = Basemap(projection='mill', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(fill_color='xkcd:light grey')

    # Iterate over each dataframe in list_of_dfs and draw scatter plot for each satellite type
    for df in list_of_dfs:
        for satellite_type, launches in selected_satellites.items():
            for launch, _ in launches.items():
                if 'L' + str(df['launch_no'][0]) in selected_satellites[satellite_type]:
                    col = launch_colour_dict[launch]
                    m.scatter(df['lons'], df['lats'], latlon=True, alpha=0.1, s=0.01, c=col)
    
    # Add parallels and meridians
    m.drawparallels(np.arange(-90., 120., 30.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 60.), labels=[0,0,0,1], fontsize=10)

    # Add legend
    patches = []
    for satellite_type, launches in selected_satellites.items():
        for launch in launches:
            color = launch_colour_dict[launch]
            patches.append(mpatches.Patch(color=color, label=satellite_type.capitalize() + ' - Launch ' + launch[1:]))
    plt.legend(handles=patches, loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    # Save the figure
    plt.savefig('output/plots/ground_tracks/g_tracks.png', bbox_inches='tight', dpi=300)
    
    # Show the plot if show is set to True
    if show:
        plt.show()

if __name__ == "__main__":
    pass
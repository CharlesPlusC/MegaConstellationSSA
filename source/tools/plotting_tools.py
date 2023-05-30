"""Set of plotting tools for the project. Most of these are based on the use of lists of pandas dataframes."""
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
import json
from typing import Dict, List, Union
from collections import defaultdict
import itertools

#local imports
from .conversions import jd_to_mjd
from .analysis_tools import TLE_arglat_dict, compute_fft, sup_gp_op_benchmark, TLE_rate_dicts

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
    Generate a plot of the power spectral density (Fourier analysis) of the differences between NORAD and operator TLEs for each launch.
    
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

def plot_diff_subplots(sats_dataframe: List[pd.DataFrame], diffs: Union[str, List[str]] = 'all', show: bool = True) -> None:
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

    legend_handles = []
    legend_labels = []

    for diff in diff_types:
        ax, xlabel, ylabel, xlims = subplot_info[diff]
        data_dict = {df['constellation'].iloc[0]: df[diff] for df in sats_dataframe_list}

        text_x_offset = 0.65  # Adjust the offset to move the text to the bottom right
        text_x = 0.95  # Start the text from the right edge
        text_y = 0.95  # Start the text from the top edge
        text_y_offset = 0.2  # Adjust the offset to spread the text in the y-direction

        for constellation, data in data_dict.items():
            # Generate histogram with automatic bins
            n, bins, patches = ax.hist(data, bins='auto', alpha=0.5, label=constellation,
                                        range=xlims, color=constellation_colour_dict[constellation.lower()])

            # Calculate statistics
            mean, std = np.mean(data), np.std(data)

            # Add statistics to the plot as text
            stats_text = "\n".join([
                f'μ {constellation}: {mean:.2f}km',
                f'σ {constellation}: {std:.2f}km'
            ])
            ax.text(text_x, text_y, stats_text, transform=ax.transAxes, horizontalalignment='right',
                    verticalalignment='top', fontsize=8)
            text_y -= text_y_offset

            if constellation not in legend_labels:
                # Add only one legend entry per unique constellation
                legend_handles.append(patches[0])  # Add the first patch from each histogram
                legend_labels.append(constellation)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, which='major')
        ax.set_yscale('log')
        ax.set_ylim(bottom=1)  # set y lower limit to 1
        ax.set_xlim(xlims)  # set x limits

    # Add a single legend with one entry per constellation
    fig.legend(legend_handles, legend_labels, loc='center right')

    plt.tight_layout(rect=(0, 0, 0.85, 1))  # Adjust the rect parameter to allocate space for the legend

    plt.savefig('output/plots/histograms/all_diffs_hist.png', dpi=300)  # Save with a single filename for all differences
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

def plot_map_diffs_smallvals_all(list_of_dfs: List[pd.DataFrame], criteria: int = 1, show: bool = False) -> None:
    """
    Plot the differences that are greater than "criteria" standard deviations from the mean for all the dataframes in list_of_dfs.
    Plot them onto a geographical map (lat/lon) and save the figure to a file.

    Parameters
    ----------
    list_of_dfs : List[pd.DataFrame]
        List of dataframes containing the data to be plotted.
    criteria : int, optional
        Differences within this many standard deviations of the mean will be retained. Defaults to 1.
    show : bool, optional
        Whether to display the plot. If set to False, the plot will only be saved. Defaults to False.

    Returns
    -------
    None

    Notes
    -----
    - This function requires the Basemap and matplotlib libraries.
    - The output plots will be saved to the 'output/plots/ground_tracks/diffs_gtrax/' directory with filenames based on constellation name and difference type.
    """
    for diff_type in diff_types:

        # group data by constellation
        grouped_dfs = defaultdict(list)
        for df in list_of_dfs:
            constellationname = df['constellation'][0]
            grouped_dfs[constellationname].append(df)

        for constellationname, df_group in grouped_dfs.items():
            # make a figure with 1 subplot
            fig, axs = plt.subplots(1, 1, figsize=(8, 5))

            # Store the differences, means, and standard deviations for each DataFrame
            diff_vals = []
            mean_diff = []
            std_diff = []
            for df in df_group:
                # calculate the mean and standard deviation of the differences for this DataFrame
                diff_vals.append(df[diff_type])
                mean_diff.append(np.mean(df[diff_type]))
                std_diff.append(np.std(df[diff_type]))

            # calculate the mean of the mean and standard deviations of the differences
            mean_mean_diff = np.mean(mean_diff)
            mean_std_diff = np.mean(std_diff)
            # find the max and min values of the differences
            max_diff = mean_mean_diff + criteria*mean_std_diff
            min_diff = mean_mean_diff - criteria*mean_std_diff

            # make a new dataframe with only the values of diff_type that are within plus or minus `criteria` standard deviations of the mean
            criteria_df_list = []
            for df in df_group:
                criteria_df_list.append(df[df[diff_type] < max_diff])
                criteria_df_list.append(df[df[diff_type] > min_diff])

            criteria_df = pd.concat(criteria_df_list, ignore_index=True)

            plt.set_cmap('seismic')
            # make the map
            m = Basemap(projection='mill', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
            m.drawcoastlines()
            m.drawcountries()
            m.drawmapboundary(fill_color='xkcd:light grey')
            
            # plot the points and colour them by diff_type, where the higher the value the darker the colour
            m.scatter(criteria_df['lons'], criteria_df['lats'], latlon=True, alpha=0.3, s=0.3, c = criteria_df[diff_type], cmap='seismic', vmin=min_diff, vmax=max_diff)

            #define mappable object
            sm = plt.cm.ScalarMappable(cmap='seismic', norm=plt.Normalize(vmin=min_diff, vmax=max_diff))
            # create colour bar
            cbar = m.colorbar(sm, location='right', pad="5%")
            # add Km label to colour bar
            cbar.set_label('Km', rotation=90, labelpad=10)
            # add parallels and meridians
            m.drawparallels(np.arange(-90., 120., 30.), labels=[1,0,0,0], fontsize=10)
            m.drawmeridians(np.arange(-180., 181., 60.), labels=[0,0,0,1], fontsize=10)

            # set the title
            if diff_type == 'cart_pos_diffs':
                axs.set_title(r'$\Delta$ 3D: ' + str(constellationname), fontsize=11)
            elif diff_type == 'h_diffs':
                axs.set_title(r'$\Delta$ H: '+ str(constellationname), fontsize=11)
            elif diff_type == 'l_diffs':
                axs.set_title(r'$\Delta$ L: '+ str(constellationname), fontsize=11)
            elif diff_type == 'c_diffs':
                axs.set_title(r'$\Delta$ C: '+ str(constellationname), fontsize=11)
            # add the mean and standard deviation of the differences to the plot
            axs.text(0.05, 0.95, 'Mean: ' + str(round(mean_mean_diff, 3)) + ' Km', transform=axs.transAxes, fontsize=10, verticalalignment='top')
            axs.text(0.05, 0.90, 'Std: ' + str(round(mean_std_diff, 3)) + ' Km', transform=axs.transAxes, fontsize=10, verticalalignment='top')

            plt.tight_layout()
                # save the figure
            plt.savefig('output/plots/ground_tracks/diffs_gtrax/' + str(constellationname) + '_' + str(diff_type) + '.png', dpi=300)

            if show == True:
                plt.show()
            plt.close()  # close the plot after saving to avoid overlapping

def plot_map_diffs_smallvals_subplot(list_of_dfs: List[pd.DataFrame], criteria: int = 1, show: bool = False) -> None:
    """
    Plot the map differences for small values using subplots.

    Parameters
    ----------
    list_of_dfs : List[pd.DataFrame]
        List of DataFrames containing the data to be plotted.
    criteria : int, optional
        Criteria for determining the maximum and minimum differences, by default 1.
    show : bool, optional
        Flag indicating whether to display the plot, by default False.

    Returns
    -------
    None

    Description
    -----------
    This function plots the map differences for small values using subplots. It takes a list of DataFrames 
    containing the data to be plotted. The differences are determined based on the specified criteria.

    The map differences are plotted as scatter plots on a map projection using the Basemap library. Each 
    subplot represents a different type of difference, and each DataFrame in the list represents a different 
    constellation. The color of the scatter points represents the difference value.

    The maximum and minimum differences are calculated based on the mean and standard deviation of the 
    differences in each constellation. Data points with differences outside the calculated range are filtered 
    out for plotting.

    The resulting plots are saved to a file and can be optionally displayed.

    Notes
    -----
    - This function requires the Basemap and matplotlib libraries.
    - The output plots are saved to the 'output/plots/ground_tracks/diffs_gtrax/all_plots.png' file.
    """
    grouped_dfs = defaultdict(list)
    for df in list_of_dfs:
        constellationname = df['constellation'][0]
        grouped_dfs[constellationname].append(df)

    num_constellations = len(grouped_dfs.keys())
    num_diff_types = len(diff_types)

    # Create figure with gridspec layout for flexible subplot arrangement
    fig = plt.figure(figsize=(8*num_diff_types, 5*num_constellations))
    gs = gridspec.GridSpec(num_constellations, num_diff_types)

    for diff_type_index, diff_type in enumerate(diff_types):
        for constellation_index, (constellationname, df_group) in enumerate(grouped_dfs.items()):
            diff_vals = []
            mean_diff = []
            std_diff = []
            for df in df_group:
                diff_vals.append(df[diff_type])
                mean_diff.append(np.mean(df[diff_type]))
                std_diff.append(np.std(df[diff_type]))

            mean_mean_diff = np.mean(mean_diff)
            mean_std_diff = np.mean(std_diff)
            max_diff = mean_mean_diff + criteria*mean_std_diff
            min_diff = mean_mean_diff - criteria*mean_std_diff

            criteria_df_list = []
            for df in df_group:
                criteria_df_list.append(df[df[diff_type] < max_diff])
                criteria_df_list.append(df[df[diff_type] > min_diff])

            criteria_df = pd.concat(criteria_df_list, ignore_index=True)

            plt.set_cmap('seismic')
            
            # Add a new subplot in the grid
            axs = fig.add_subplot(gs[constellation_index, diff_type_index])
            
            m = Basemap(projection='mill', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=axs)
            m.drawcoastlines()
            m.drawcountries()
            m.drawmapboundary(fill_color='xkcd:light grey')
            m.scatter(criteria_df['lons'], criteria_df['lats'], latlon=True, alpha=0.3, s=0.3, c = criteria_df[diff_type], cmap='seismic', vmin=min_diff, vmax=max_diff)

            sm = plt.cm.ScalarMappable(cmap='seismic', norm=plt.Normalize(vmin=min_diff, vmax=max_diff))
            cbar = m.colorbar(sm, location='right', pad="5%")
            cbar.set_label('Km', rotation=90, labelpad=5)
            m.drawparallels(np.arange(-90., 120., 30.), labels=[1,0,0,0], fontsize=14)
            m.drawmeridians(np.arange(-180., 181., 60.), labels=[0,0,0,1], fontsize=14)

            # set the title
            if diff_type == 'cart_pos_diffs':
                axs.set_title(r'$\Delta$ 3D: ' + str(constellationname), fontsize=14)
            elif diff_type == 'h_diffs':
                axs.set_title(r'$\Delta$ H: '+ str(constellationname), fontsize=14)
            elif diff_type == 'l_diffs':
                axs.set_title(r'$\Delta$ L: '+ str(constellationname), fontsize=14)
            elif diff_type == 'c_diffs':
                axs.set_title(r'$\Delta$ C: '+ str(constellationname), fontsize=14)

            axs.text(0.05, 0.95, 'Mean: ' + str(round(mean_mean_diff, 3)) + ' Km', transform=axs.transAxes, fontsize=14, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
            axs.text(0.05, 0.90, 'Std: ' + str(round(mean_std_diff, 3)) + ' Km', transform=axs.transAxes, fontsize=14, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

    plt.tight_layout()
    plt.savefig('output/plots/ground_tracks/diffs_gtrax/all_plots.png', dpi=300)
    if show == True:
        plt.show()
    plt.close()

def benchmark_plot() -> None:
    """
    Plot the benchmarking analysis results.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Description
    -----------
    This function plots differences in orbital elements H, C, L, and 3D positions between 
    NORAD (GP), supplementary (SUP), and operator (Op) TLE data for the satellites whose data is in "external/ephem_TLE_compare". 

    The differences are plotted as a function of time in Modified Julian Date (MJD) format. 
    Vertical dotted lines represent the epochs of the respective TLE data. 

    Each NORAD ID is plotted in a separate subplot column with four rows of plots for H, C, L, and 3D differences.

    Notes
    -----
    This function expects that the required data has been precomputed by the function `sup_gp_op_benchmark()`.
    """
    all_triple_ephems, all_sup_tle_epochs, all_gp_tle_epochs, gp_list = sup_gp_op_benchmark()

    #convert all_sup_tle_epochs and all_gp_tle_epochs to mjd_time
    for i in range(len(all_sup_tle_epochs)):
        all_sup_tle_epochs[i] = [i - 2400000.5 for i in all_sup_tle_epochs[i]]

    for i in range(len(all_gp_tle_epochs)):
        all_gp_tle_epochs[i] = [i - 2400000.5 for i in all_gp_tle_epochs[i]]

    for i in range(len(all_triple_ephems)):
        print(len(all_triple_ephems[i]))
        #assert they are all the same length
        # assert len(all_triple_ephems[i]) == len(all_sup_tle_epochs[i]) == len(all_gp_tle_epochs[i])

    #get a list of NORAD IDs from the GP ephemeris files
    NORAD_IDs = []
    for gp_path in gp_list:
        NORAD_IDs.append(gp_path[-9:-4])

    ### Plot all the data in all_triple_ephems (4 subplot plots per NORAD dataset so 4*len(all_triple_ephems) plots total)
    fig, ax = plt.subplots(4, len(all_triple_ephems), figsize=(4*len(all_triple_ephems), 6))

    GP_colour = 'xkcd:blue'
    SUP_colour = 'xkcd:coral'

    fonts = 12
    for i in range(len(all_triple_ephems)):
        # for each NORAD ID plot h, c, l, and 3D position differences against jd_time for both GP and SUP
        #only include y-axis labels for the first plot
        if i == 0:
            ax[0, i].set_ylabel(r'$\Delta$ H (km)', fontsize=fonts)
        # else remove the tick labels
        else:
            ax[0, i].set_yticklabels([])
        ax[0, i].plot(all_triple_ephems[i]['mjd_time'].values, all_triple_ephems[i]['gp_h_diff'].values, label='GP', color=GP_colour)
        ax[0, i].plot(all_triple_ephems[i]['mjd_time'].values, all_triple_ephems[i]['sup_h_diff'].values, label='SUP', color=SUP_colour)
        # add vertical dotted thin lines for each element in the list all_sup_tle_epochs[i]
        for j in range(len(all_sup_tle_epochs[i])):
            ax[0, i].axvline(x=all_sup_tle_epochs[i][j], color=SUP_colour, linestyle='dotted', linewidth=2)
        # same for GP
        for j in range(len(all_gp_tle_epochs[i])):
            ax[0, i].axvline(x=all_gp_tle_epochs[i][j], color=GP_colour, linestyle='dotted', linewidth=2)
        ax[0, i].grid(True)
        #remove x-axis labels
        ax[0, i].set_xticklabels([])
        ax[0, i].set_title('NORAD ID: ' + NORAD_IDs[i], fontsize=fonts)
        # for the y labels to be from -0.25 to 0.25
        ax[0, i].set_ylim(-0.5, 0.5)
        # set x-axis limits to be the max and min of the mjd_time column
        ax[0, i].set_xlim(min(all_triple_ephems[i]['mjd_time']), max(all_triple_ephems[i]['mjd_time']))

        if i == 0:
            ax[1, i].set_ylabel(r'$\Delta$ C (km)', fontsize=fonts)
        # else remove the tick labels
        else:
            ax[1, i].set_yticklabels([])
        ax[1, i].plot(all_triple_ephems[i]['mjd_time'].values, all_triple_ephems[i]['gp_c_diff'].values, label='GP', color=GP_colour)
        ax[1, i].plot(all_triple_ephems[i]['mjd_time'].values, all_triple_ephems[i]['sup_c_diff'].values, label='SUP', color=SUP_colour)
        # add vertical dotted thin lines for each element in the list all_sup_tle_epochs[i]
        for j in range(len(all_sup_tle_epochs[i])):
            ax[1, i].axvline(x=all_sup_tle_epochs[i][j], color=SUP_colour, linestyle='dotted', linewidth=2)
        # same for GP
        for j in range(len(all_gp_tle_epochs[i])):
            ax[1, i].axvline(x=all_gp_tle_epochs[i][j], color=GP_colour, linestyle='dotted', linewidth=2)
        ax[1, i].set_xticklabels([])
        ax[1, i].grid(True)
        # set x-axis limits to be the max and min of the mjd_time column
        ax[1, i].set_xlim(min(all_triple_ephems[i]['mjd_time']), max(all_triple_ephems[i]['mjd_time']))
        # for the y labels to be from -0.6 to 0.6
        # ax[1, i].set_ylim(-0.5, 0.5)

        if i == 0:
            ax[2, i].set_ylabel(r'$\Delta$ L (km)', fontsize=fonts)
        # else remove the tick labels
        else:
            ax[2, i].set_yticklabels([])
        ax[2, i].plot(all_triple_ephems[i]['mjd_time'].values, all_triple_ephems[i]['gp_l_diff'].values, label='GP', color=GP_colour)
        ax[2, i].plot(all_triple_ephems[i]['mjd_time'].values, all_triple_ephems[i]['sup_l_diff'].values, label='SUP', color=SUP_colour)
        # add vertical dotted thin lines for each element in the list all_sup_tle_epochs[i]
        for j in range(len(all_sup_tle_epochs[i])):
            ax[2, i].axvline(x=all_sup_tle_epochs[i][j], color=SUP_colour, linestyle='dotted', linewidth=2)
        # same for GP
        for j in range(len(all_gp_tle_epochs[i])):
            ax[2, i].axvline(x=all_gp_tle_epochs[i][j], color=GP_colour, linestyle='dotted', linewidth=2)
        ax[2, i].set_xticklabels([])
        ax[2, i].grid(True)
        ax[2, i].set_ylim(-10, 10)
        # set x-axis limits to be the max and min of the mjd_time column
        ax[2, i].set_xlim(min(all_triple_ephems[i]['mjd_time']), max(all_triple_ephems[i]['mjd_time']))

        if i == 0:
            ax[3, i].set_ylabel(r'$\Delta$ 3D (km)', fontsize=fonts)
        # else remove the tick labels
        else:
            ax[3, i].set_yticklabels([])
        ax[3, i].plot(all_triple_ephems[i]['mjd_time'].values, all_triple_ephems[i]['gp_cart_pos_diff'].values, label='GP', color=GP_colour)
        ax[3, i].plot(all_triple_ephems[i]['mjd_time'].values, all_triple_ephems[i]['sup_cart_pos_diff'].values, label='SUP', color=SUP_colour)
        # add vertical dotted thin lines for each element in the list all_sup_tle_epochs[i]
        for j in range(len(all_sup_tle_epochs[i])):
            ax[3, i].axvline(x=all_sup_tle_epochs[i][j], color=SUP_colour, linestyle='dotted', linewidth=2)
        # same for GP
        for j in range(len(all_gp_tle_epochs[i])):
            ax[3, i].axvline(x=all_gp_tle_epochs[i][j], color=GP_colour, linestyle='dotted', linewidth=2)
        # calculate the RMS of the differences
        gp_rms = np.sqrt(np.mean(np.square(all_triple_ephems[i]['gp_cart_pos_diff'])))
        sup_rms = np.sqrt(np.mean(np.square(all_triple_ephems[i]['sup_cart_pos_diff'])))
        gp_mean = np.mean(all_triple_ephems[i]['gp_cart_pos_diff'])
        sup_mean = np.mean(all_triple_ephems[i]['sup_cart_pos_diff'])
        # add a little text box with the RMS values
        textstr = 'GP mean: ' + str(round(gp_mean, 2)) + ' km \n SUP mean: ' + str(round(sup_mean, 2)) + ' km'
        # props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        # ax[3, i].text(0.45, 0.9, textstr, transform=ax[3, i].transAxes, fontsize=fonts, verticalalignment='top', bbox=props)
        ax[3, i].set_title(textstr, fontsize=fonts)
        ax[3, i].grid(True)
        ax[3, i].set_ylim(0, 8)
        # set ticks every 2.5 on the y axis
        ax[3, i].yaxis.set_ticks(np.arange(0, 8.1, 2))
        # set x-axis limits to be the max and min of the mjd_time column
        ax[3, i].set_xlim(min(all_triple_ephems[i]['mjd_time']), max(all_triple_ephems[i]['mjd_time']))
        # include a common x label for the bottom row make it not scientific notation
        ax[3, i].set_xlabel('MJD Time', fontsize=fonts)
        ax[3, i].ticklabel_format(style='plain', axis='x', scilimits=(0,0)) 

    # include one legend for the whole figure but only for the lines in the first plot
    handles, labels = ax[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', ncol=1, fontsize=fonts)

    plt.tight_layout()
    plt.savefig('output/plots/benchmark/SUPvsGPvsOp_06022023.png', dpi=300)

    plt.show()

def plot_arglat_analysis(show: bool = False) -> None:
    """
    Plot argument of latitude analysis of TLEs.

    This function fetches the TLE data (converted to argument of latitude values) 
    from two sources for OneWeb and Starlink satellites and generates a set of histograms 
    to visualize the distribution of the argument of latitude for these satellites. 
    It saves the generated plot to an external file and optionally displays the plot.

    Parameters
    ----------
    show : bool, optional
        If True, the generated plot is displayed. The default is False.

    Returns
    -------
    None
    """
    GP_arglats = TLE_arglat_dict(selected_satellites='external/selected_satellites.json', tle_folder = 'external/NORAD_TLEs/')
    SUP_arglats = TLE_arglat_dict(selected_satellites='external/selected_satellites.json', tle_folder = 'external/SUP_TLEs/')

    # Flatten list of lists
    oneweb_gp = [value for key, value in GP_arglats["oneweb"].items()]
    oneweb_gp_flat = [arglat for sublist in oneweb_gp for arglat in sublist]
    
    starlink_gp = [value for key, value in GP_arglats["starlink"].items()]
    starlink_gp_flat = [arglat for sublist in starlink_gp for arglat in sublist]
    
    oneweb_sup = [value for key, value in SUP_arglats["oneweb"].items()]
    oneweb_sup_flat = [arglat for sublist in oneweb_sup for arglat in sublist]
    
    starlink_sup = [value for key, value in SUP_arglats["starlink"].items()]
    starlink_sup_flat = [arglat for sublist in starlink_sup for arglat in sublist]

    # Calculate the total number of data points in each list
    ow_norad_total = len(oneweb_gp_flat)
    ow_sup_total = len(oneweb_sup_flat)
    sl_norad_total = len(starlink_gp_flat)
    sl_sup_total = len(starlink_sup_flat)

    # Create weights for each dataset to normalize the histogram counts
    ow_norad_weights = np.ones_like(oneweb_gp_flat) / ow_norad_total
    ow_sup_weights = np.ones_like(oneweb_sup_flat) / ow_sup_total
    sl_norad_weights = np.ones_like(starlink_gp_flat) / sl_norad_total
    sl_sup_weights = np.ones_like(starlink_sup_flat) / sl_sup_total

    fig, axs = plt.subplots(2, 2, figsize=(8, 5), sharex='col', sharey='none', 
                            gridspec_kw={'height_ratios': [1, 3], 'wspace': 0.05, 'hspace': 0.15})

    axs[0, 0].hist(oneweb_gp_flat, bins=30, weights=ow_norad_weights, alpha=0.5, label='NORAD TLEs', color="xkcd:azure")
    axs[0, 0].hist(oneweb_sup_flat, bins=30, weights=ow_sup_weights, alpha=0.5, label='SUP TLEs', color="xkcd:blue")
    axs[0, 0].set_ylim(0.7, 0.8)

    axs[1, 0].hist(oneweb_gp_flat, bins=30, weights=ow_norad_weights, alpha=0.5, label='NORAD TLEs', color="xkcd:azure")
    axs[1, 0].hist(oneweb_sup_flat, bins=30, weights=ow_sup_weights, alpha=0.5, label='SUP TLEs', color="xkcd:blue")
    axs[1, 0].set_ylim(0, 0.3)

    axs[0, 1].hist(starlink_gp_flat, bins=30, weights=sl_norad_weights, alpha=0.5, label='NORAD TLEs', color="xkcd:orange")
    axs[0, 1].hist(starlink_sup_flat, bins=30, weights=sl_sup_weights, alpha=0.5, label='SUP TLEs', color="xkcd:coral")
    axs[0, 1].set_ylim(0.7, 0.8)

    axs[1, 1].hist(starlink_gp_flat, bins=30, weights=sl_norad_weights, alpha=0.5, label='NORAD TLEs', color="xkcd:orange")
    axs[1, 1].hist(starlink_sup_flat, bins=30, weights=sl_sup_weights, alpha=0.5, label='SUP TLEs', color="xkcd:coral")
    axs[1, 1].set_ylim(0, 0.3)

        # Additional plot styling
    for ax in axs.flat:
        ax.set_xlim(-1, 361)
        ax.set_xticks(np.arange(0, 361, 60))
        # make the ticks diagonal so they don't overlap
        ax.tick_params(axis='x', rotation=35)
        ax.grid()
        
    for ax in axs[1,:]:
        ax.set_xlabel('Argument of Latitude (deg)')
        ax.set_ylabel('')
        
    for ax in axs[:,0]:
        ax.set_ylabel('Proportion of TLEs')
    
    #top left subplot should have a blank y label
    axs[0,0].set_ylabel('')
        
    axs[0,0].legend(loc='upper right', title='OneWeb')
    axs[0,1].legend(loc='upper right', title='Starlink')

    # # Removing xticklabels for the upper plots
    # for ax in axs[0, :]:
    #     ax.set_xticklabels([])

    # Removing yticklabels for the second column plots
    for ax in axs[:, 1]:
        ax.set_yticklabels([])

    plt.tight_layout()
    plt.savefig('output/plots/TLE_production/tle_arglat_hist.png', dpi=300, bbox_inches='tight')
    if show: 
        plt.show()

def plot_tle_rate_analysis(show: bool = False) -> None:
    """
    Plot TLE latency (/production rate) analysis.

    This function fetches the TLE data for OneWeb and Starlink satellites 
    and generates a set of histograms to visualize the latency between successive TLEs 
    for the satellites in these constellations. It saves the generated plot to an external 
    file and optionally displays the plot.

    Parameters
    ----------
    show : bool, optional
        If True, the generated plot is displayed. The default is False.

    Returns
    -------
    None
    """
    NORAD_rate_dicts = TLE_rate_dicts(tle_folder = 'external/NORAD_TLEs/')
    SUP_rate_dicts = TLE_rate_dicts(tle_folder = 'external/SUP_TLEs/')

    NORAD_ow = NORAD_rate_dicts[0] #Oneweb
    NORAD_sl = NORAD_rate_dicts[1] #Starlink

    SUP_ow = SUP_rate_dicts[0] #Oneweb
    SUP_sl = SUP_rate_dicts[1] #Starlink

    # Flatten the values of the four dictionaries into lists
    #TODO: there is definitely a better way of doing this
    NORAD_ow_flat = list(itertools.chain.from_iterable(NORAD_ow.values()))
    NORAD_sl_flat = list(itertools.chain.from_iterable(NORAD_sl.values()))
    SUP_ow_flat = list(itertools.chain.from_iterable(SUP_ow.values()))
    SUP_sl_flat = list(itertools.chain.from_iterable(SUP_sl.values()))

    # now replicate the two plots above but put them side by side in subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,5))

    bins = np.arange(0, 40, 2)

    ax1.hist(NORAD_ow_flat, bins=bins, alpha=0.5, label='NORAD TLEs', color="xkcd:azure")
    ax1.hist(SUP_ow_flat, bins=bins, alpha=0.5, label='SUP TLEs', color="xkcd:blue")
    ax1.set_xlim(0,30)
    ax1.set_ylim(0,25000)
    #add ticks every 2 hours
    ax1.set_xticks(np.arange(0, 30.1, 2))
    # rotate the ticks so they don't overlap
    ax1.tick_params(axis='x', rotation=35)
    
    # calculate the mean and standard deviation of the NORAD TLEs and the SUP TLEs and add them to the plot on the middle of the right 
    # side of the plot
    # ax1.text(2, 20000, 'NORAD TLEs mean: {:.2f} hours '.format(np.mean(NORAD_ow_flat)))
    # ax1.text(2, 18000, 'NORAD TLEs std: {:.2f} hours '.format(np.std(NORAD_ow_flat)))
    # ax1.text(2, 16000, 'SUP TLEs mean: {:.2f} hours '.format(np.mean(SUP_ow_flat)))
    # ax1.text(2, 14000, 'SUP TLEs std: {:.2f} hours '.format(np.std(SUP_ow_flat)))
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel('TLE Latency (Hours/TLE)')
    ax1.set_ylabel('Number of TLEs')
    ax1.set_title('OneWeb')

    bins = np.arange(0, 40, 2)
    ax2.hist(NORAD_sl_flat, bins=bins, alpha=0.5, label='NORAD TLEs', color="xkcd:orange")
    ax2.hist(SUP_sl_flat, bins=bins, alpha=0.5, label='SUP TLEs', color="xkcd:coral")
    ax2.set_xlim(0,30)
    ax2.set_ylim(0,25000)
    #add ticks every 2 hours
    ax2.set_xticks(np.arange(0, 30.1, 2))
    # rotate the ticks so they don't overlap
    ax2.tick_params(axis='x', rotation=35)
    #remove the y ticks
    ax2.set_yticklabels([])

    # calculate the mean and standard deviation of the NORAD TLEs and the SUP TLEs and add them to the plot on the middle of the right
    # side of the plot
    # ax2.text(2, 20000, 'NORAD TLEs mean: {:.2f} hours '.format(np.mean(NORAD_sl_flat)))
    # ax2.text(2, 18000, 'NORAD TLEs std: {:.2f} hours '.format(np.std(NORAD_sl_flat)))
    # ax2.text(2, 16000, 'SUP TLEs mean: {:.2f} hours '.format(np.mean(SUP_sl_flat)))
    # ax2.text(2, 14000, 'SUP TLEs std: {:.2f} hours '.format(np.std(SUP_sl_flat)))

    ax2.legend()
    ax2.grid()
    ax2.set_xlabel('TLE Latency (Hours/TLE)')
    ax2.set_ylabel('')
    ax2.set_title('Starlink')

    plt.tight_layout()

    # save the figure
    fig.savefig('output/plots/TLE_production/TLE_Latencies.png', dpi=300)

    if show: 
        plt.show()

if __name__ == "__main__":
    pass
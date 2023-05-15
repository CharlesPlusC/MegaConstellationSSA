"""Set of plotting tools for the project. Most of these are based on the use of lists of pandas dataframes."""

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
import json

mpl.rcParams['font.size'] = 11

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
import json

mpl.rcParams['font.size'] = 11

def plot_altitude_timeseries(dfs, json_filepath='external/satellite_lists.json', show=False):
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
                            if launch in ['L4', 'L28']:
                                colours.append('xkcd:blue')
                            elif launch in ['L5', 'L36']:
                                colours.append('xkcd:azure')
                            elif launch in ['L6', 'L30']:
                                colours.append('xkcd:light blue')

    fig, ax = plt.subplots(figsize=(8, 5))
    for i in range(len(eph_alts_sup)):
        ax.scatter(x=times[i], y=eph_alts_sup[i], color=colours[i], label=launch_nos[i], alpha=0.2, s=0.1)
    ax.set_title(constellation)
    ax.set_ylabel('Altitude (Km)')
    ax.set_xlabel('Time (MJD)')
    ax.set_xlim(59390, 59940) # window of interest for this project
    ax.legend(loc='upper right')
    ax.grid()
    plt.tight_layout()
    plt.savefig('altitude_tseries_all.png', dpi=400)
    if show:
        plt.show()

import os
import sys
from tools.analysis_tools import NORAD_vs_SUP_TLE_analysis
from tools.tletools import NORAD_TLE_History, load_satellite_lists

satellite_lists = load_satellite_lists()

# Access the lists using the satellite_lists dictionary
norad_OW_L5 = satellite_lists["oneweb"]["L5"]
norad_OW_L4 = satellite_lists["oneweb"]["L4"]
norad_OW_L6 = satellite_lists["oneweb"]["L6"]

norad_SL_L28 = satellite_lists["starlink"]["L28"]
norad_SL_L30 = satellite_lists["starlink"]["L30"]
norad_SL_L36 = satellite_lists["starlink"]["L36"]

#merge the lists above into one list
SL_norads = norad_SL_L28 + norad_SL_L36 + norad_SL_L30
OW_norads = norad_OW_L4 + norad_OW_L5 + norad_OW_L6

#convert the list of ints to strings for the function
SL_norads = [str(i) for i in SL_norads]
OW_norads = [str(i) for i in OW_norads]

# fetch the available TLE archive for the specified NORAD IDs
NORAD_TLE_History(SL_norads, 'starlink') 
NORAD_TLE_History(OW_norads, 'oneweb')

# NORAD_vs_SUP_TLE_analysis(NORADS = SL_norads) #run the function on the SL satellites
# NORAD_vs_SUP_TLE_analysis(NORADS = OW_norads) #run the function on the OW satellites

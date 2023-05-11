from source.tools.analysis_tools import NORAD_vs_SUP_TLE_analysis

#Manually define the selected satellites for reproducability
norad_OW_L5 = [48042, 48043, 48044, 48045, 48046, 48047, 48048, 48049, 48050, 48051]
norad_OW_L4 = [47269,47280,47258,47279,47277,47273,47260,47264,47266,47288]
norad_OW_L6 = [48210, 48211, 48212, 48213, 48214, 48215, 48216, 48217, 48218, 48219]

norad_SL_L28 = [48430, 48434,48432,48431,48430,48435,48436,48439,48437,48433]
norad_SL_L30 =  [48638,48639,48640,48641,48642,48643,48644,48645,48646,48647]
norad_SL_L36 = [50803, 50804, 50805, 50806, 50807, 50808, 50809, 50810, 50811, 50812]

#merge the lists above into one list
SL_norads = norad_SL_L28 + norad_SL_L36 + norad_SL_L30
OW_norads = norad_OW_L4 + norad_OW_L5 + norad_OW_L6

#convert the list of ints to strings for the function
SL_norads = [str(i) for i in SL_norads]
OW_norads = [str(i) for i in OW_norads]

NORAD_vs_SUP_TLE_analysis(NORADS = SL_norads) #run the function on the SL satellites
NORAD_vs_SUP_TLE_analysis(NORADS = OW_norads) #run the function on the OW satellites
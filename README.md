# MegaConstellationSSA
Code repository accompanying the paper: On the Limits of Current Practices in Uncooperative Space Surveillance: Analysis of Mega-Constellation Data Time-Series

## Instructions for use
1. Clone the repository
2. Install the requirements using the megeaconst_env.yml file:
``` bash
(conda env create -f megaconst_env.yml)
```
3. Activate the environment:
``` bash
 (conda activate megaconst)
```
4. Create a file called SLTrack.ini and put it in the root directory of the repository. The file should contain your Space-track.org username and password in the following format:
``` bash
[configuration]
username = your_email@email.com
password = YourPassword
```
5. Ensure that this file is in the .gitignore file so that it is not uploaded to the repository (it should be by default)
6. Run the code in main.py for the full analysis

## Data
The data for the SUP_TLEs and NORAD_TLEs are already provided to save time. 
If you wish to download the data yourself, you can do so for the NORAD TLEs by running the getdata.py script. Note that there is currently no way of downloading SUP_TLE data programatically so you will have to do this manually from celestrak.org.

## Issues
If you have any issues with the code, please raise an issue on this repository and I will try to get back to you as soon as possible.

## License: GNU GPLv3

## Envs
to save your current enviro
``` bash
conda env export > megaconst_env.yml
```
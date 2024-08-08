# Library Imports
import numpy as np


###### Unit Vectors ######
i_index = 0
j_index = 1
k_index = 2

i_unit = np.array([1, 0, 0])
j_unit = np.array([0, 1, 0])                 
k_unit = np.array([0, 0, 1])


###### Unit constants ######
GRAVITATIONAL_PARAMETER = 1  

TU_TO_SEC = (1 / 1.239446309e-3)

KM_TO_CANONICAL = 1 / 6371
SEC_TO_TU = 806.8118744
FT_TO_CANONICAL = 1 / 20902302

EARTH_ECCENTRICITY = 0.08182
EARTH_ANGULAR_VELOCITY = np.array([0, 0, 0.0588336565]) 


###### Test Case ######
LAUNCH_SITE = { "latitude": 16.45,
                "longitude": 169.32,
                "altitude": 5 }

RADAR_TRACKING_SITE = { "latitude": 52.45,
                        "longitude": 174.05,
                        "altitude": 52,
                        "day": 104,
                        "universal_time": "0600:00" }

RADAR_MEASUREMENTS = { "range": 186.613,
                       "range_p": 4.012,
                       "elevation": 56.95,
                       "elevation_p": -1.92,
                       "azimuth": 152.44,
                       "azimuth_p": 1.09 }

# test case settings
FIRST_REC_T = 10
REC_T_INCREMENT_SIZE = 5
NUM_REC_TIMES = 40
FIRST_TOF = 5
TOF_INCREMENT_SIZE = 5
NUM_TOFS = 14 

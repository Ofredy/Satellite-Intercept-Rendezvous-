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

FT_TO_CANONICAL = 1 / 20902302

EARTH_ECCENTRICITY = 0.08182
EARTH_ANGULAR_VELOCITY = np.array([0, 0, 0.0588336565]) 

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


###### Testcases ######
test_case_positions = [[0, 1, 0],
                       [0, 0, -0.5],
                       [0.3, 1, 0],
                       [0.5, 0.7, 0.8],
                       [0.025917, -0.150689, 1.138878],
                       [-0.5, 0, 0]]

test_case_velocities = [[ 0, 0, 1],
                        [ 0, 2, 0],
                        [ 3, 0, 0],
                        [ 0, 0.1, 0.9],
                        [ 0.000361, 0.001074, 0.002177],
                        [ 0, 1.999, 0]]

test_case_t = [ np.pi,
                10**6,
                5,
                -20,
                1.5,
                10**3 ]
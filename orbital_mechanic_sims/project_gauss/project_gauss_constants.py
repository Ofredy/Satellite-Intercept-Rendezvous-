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
test_case_r1 = [ [0.5, 0.6, 0.7],
                 [0.3, 0.7, 0.4],
                 [0.5, 0.6, 0.7],
                 [-0.2, 0.6, 0.3],
                 [1, 0, 0],
                 [-0.4, 0.6, -1.201]]

test_case_r2 = [ [0, -1, 0],
                 [0.6, -1.4, 0.8],
                 [0, 1, 0],
                 [0.4, 1.2, 0.6],
                 [0, 1, 0],
                 [0.2, -0.3, 0.6]]

test_case_t = [ 20,
                5,
                1.2,
                50,
                0.0001,
                5 ]

test_case_dm = [ -1,
                 1,
                 -1,
                 1,
                 1,
                 1 ]
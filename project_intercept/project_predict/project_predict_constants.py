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
TU_TO_SEC = (1 / 1.239446309e-3)


###### Testcases ######
test_case_positions = [ [-.1, 1, 0],
                        [0, 0, 2],
                        [.49, .48, .9],
                        [0, 4, 0],
                        [0, 0, 2],
                        [-2.414, -2.414, 0],
                        [0, 2.1, .001],
                        [0, 0, 530],
                        [-65.62, 22.9, 0]]

test_case_velocities = [ [-1.2, -.01, 0],
                         [ 0, -.49, -.1],
                         [ 0, 0, 1.01],
                         [ -.5, -.5, 0],
                         [ .8, 0, .6],
                         [ .707, .293, 0],
                         [ -.703, -.703, .001],
                         [ -.00001, -.05, -1],
                         [ .01745, .000305, 0]]
####################### File Description #######################
#     File used for Project Predict.
#
#     For a number of unidentified space objects the three components of
#     vector position and velocity are generated from radar observations (see
#     Project SITE/TRACK). From this information you are required to find:
#
#       a. The type of trajectory (circular, rectilinear, elliptical, parabolic,
#       hyperbolic).
#
#       b. Position and velocity vectors at impact or closest approach in the
#       geocentric-equatorial coordinate system.
#
#       c. Time for object to go from its observed position to point of
#       impact or closest approach.
#
#       d. The total change in true anomaly from the observed position to
#       impact or point of closest approach. Check to see if trajectory is going
#       away from the earth. 


# Library Inputs
import numpy as np


class ProjectPredict:

    def __init__(self, position_vector: np.array, velocity_vector: np.array):
        
        self.position_vector = position_vector
        self.velocity_vector = velocity_vector

        self._find_orbital_parameters()


    def _find_orbital_parameters(self):

        pass
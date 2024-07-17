####################### File Description #######################
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


# System Imports
import math

# Library Imports
import numpy as np

# Our Imports
from project_predict_constants import *


class Satellite:

    def __init__(self, position_vector: list, velocity_vector: list):
        
        self.r = np.expand_dims(np.array(position_vector), axis=0)
        self.v = np.expand_dims(np.array(velocity_vector), axis=0)

        self._find_orbital_parameters()

        self._find_trajectory_type()

        self._find_impact_or_closest_approach()

    def _find_orbital_parameters(self):

        # angular momentum
        self.h = np.cross(self.r, self.v)

        # semi-latus rectus 
        self.p = (np.linalg.norm(self.h)**2) / GRAVITATIONAL_PARAMETER

        # eccentricity
        self.e = (1/GRAVITATIONAL_PARAMETER) * ((np.linalg.norm(self.v)**2 - (GRAVITATIONAL_PARAMETER/np.linalg.norm(self.r))) * self.r - (np.vdot(self.r, self.v)) * self.v)

        # inclination
        self.i = math.acos(np.vdot(self.h, k_unit) / np.linalg.norm(self.h))

        # node vector
        self.n = np.cross(k_unit, self.h)

        # longitude of ascending node
        self.omega = math.acos(np.vdot(self.n, i_unit) / np.linalg.norm(self.n))

        if self.n[0][j_index] < 0:
            self.omega = 2 * math.pi - self.omega

        # argument of periapsis
        self.w = math.acos(np.vdot(self.n, self.e) / (np.linalg.norm(self.n) * np.linalg.norm(self.e)))

        if self.e[0][k_index] < 0:
            self.w = 2 * math.pi - self.w

        # true anomaly
        self.v_o = math.acos(np.vdot(self.e, self.r) / (np.linalg.norm(self.e) * np.linalg.norm(self.r)))

        if np.vdot(self.r, self.v) < 0:
            self.v_0 = 2 * math.pi - self.v_o

    def _find_trajectory_type(self):

        e_norm = np.linalg.norm(self.e)

        if e_norm == 0:
            self.trajectory_type = "circular"

        elif 0 < e_norm and e_norm < 1:
            self.trajectory_type = "elliptical"
            
        elif e_norm == 1:
            self.trajectory_type = "parabolic"

        elif e_norm > 1:
            self.trajectory_type = "hyperbolic"

        elif np.linalg.norm(self.h) == 0:
            self.trajectory_type = "rectilinear"

        else:
            self.trajectory_type = "undefined"

    def _find_impact_or_closest_approach(self):

        pass


if __name__ == "__main__":

    satellite = Satellite( [3*math.sqrt(3)/4, 3/4, 0], [-1/(2*math.sqrt(2)), math.sqrt(3)/(2*math.sqrt(2)), 1/math.sqrt(2)] )

    
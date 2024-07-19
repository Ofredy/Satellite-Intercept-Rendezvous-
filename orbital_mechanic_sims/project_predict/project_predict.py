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
        
        self.r = np.array(position_vector)
        self.v = np.array(velocity_vector)

        self._find_orbital_parameters()

        self._find_trajectory_type()

        #self._find_impact_or_closest_approach()

    def _find_orbital_parameters(self):

        # angular momentum
        self.h = np.cross(self.r.flatten(), self.v.flatten()).reshape(3, 1)

        # semi-latus rectus 
        self.p = (np.linalg.norm(self.h)**2) / GRAVITATIONAL_PARAMETER

        # eccentricity
        self.e = (1/GRAVITATIONAL_PARAMETER) * ((np.linalg.norm(self.v)**2 - (GRAVITATIONAL_PARAMETER/np.linalg.norm(self.r))) * self.r - (np.vdot(self.r, self.v)) * self.v)

        # inclination
        self.i = math.acos(np.vdot(self.h, k_unit) / np.linalg.norm(self.h))

        # node vector
        self.n = np.cross(k_unit.flatten(), self.h.flatten()).reshape(3, 1)

        # longitude of ascending node
        self.omega = math.acos(np.vdot(self.n, i_unit) / np.linalg.norm(self.n))

        if self.n[j_index][0] < 0:
            self.omega = 2 * math.pi - self.omega

        # argument of periapsis
        self.w = math.acos(np.vdot(self.n, self.e) / (np.linalg.norm(self.n) * np.linalg.norm(self.e)))

        if self.e[k_index][0] < 0:
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

        # solving for position vector at perigee
        true_anomaly_perigee = 0
        perigee_radius = self.p / ( 1 + self.e * math.cos(true_anomaly_perigee))

        r_p_perifocal_coordinate = np.expand_dims(np.array([perigee_radius*math.cos(true_anomaly_perigee), perigee_radius*math.sin(true_anomaly_perigee), 0]))
        self._solve_for_perifocal_to_ijk_matrix()

        r_p = self.perifocal_to_ijk_matrix * r_p_perifocal_coordinate

        # checking to see if an impact occured
        if np.linalg.norm(r_p) < 1:
            self.impact = True
            self._solve_for_impact_point()

    def _solve_for_perifocal_to_ijk_matrix(self):
        
        self.perifocal_to_ijk_matrix = np.array([[ math.cos(self.omega)*math.cos(self.w)-math.sin(self.omega)*math.sin(self.w)*math.cos(self.i), -1*math.cos(self.omega)]])

    def _solve_for_impact_point(self):

        pass

if __name__ == "__main__":

    satellite = Satellite( [[3*math.sqrt(3)/4], [3/4], [0]], [[-1/(2*math.sqrt(2))], [math.sqrt(3)/(2*math.sqrt(2))], [1/math.sqrt(2)]] )

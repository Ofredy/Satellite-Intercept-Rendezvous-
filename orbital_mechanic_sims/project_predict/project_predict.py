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
        self._find_perigee_or_impact()
        self._find_time_of_flight()

    def _find_orbital_parameters(self):

        # angular momentum
        self.h = np.cross(self.r, self.v)

        # semi-latus rectus 
        self.p = (np.linalg.norm(self.h)**2) / GRAVITATIONAL_PARAMETER

        # eccentricity
        self.e = (np.linalg.norm(self.v)**2 - (GRAVITATIONAL_PARAMETER/np.linalg.norm(self.r))) * self.r - np.dot(self.r, self.v) * self.v

        # semi-major axis
        self.a = self.p / ( 1 - np.linalg.norm(self.e)**2 )

        # inclination
        self.i = np.arccos(np.dot(self.h, k_unit) / np.linalg.norm(self.h))

        # node vector
        self.n = np.cross(k_unit, self.h)

        # longitude of ascending node and argument of periapsis
        if abs(self.i) < 1e-10 or abs(self.i - math.pi) < 1e-10:
            # equatorial orbit -> omega & w become undefined so we set them to zero for future calculations
            self.omega = 0
            self.w = 0

        else:
            # non equatorial orbit
            self.omega = np.arccos(np.dot(self.n, i_unit) / np.linalg.norm(self.n))

            if self.n[j_index] < 0:
                self.omega = 2 * math.pi - self.omega

            self.w = np.arccos(np.dot(self.n, self.e) / (np.linalg.norm(self.n) * np.linalg.norm(self.e)))

            if self.e[k_index] < 0:
                self.w = 2 * math.pi - self.w

        # true anomaly
        self.nu_0 = np.arccos(np.dot(self.e, self.r) / (np.linalg.norm(self.e) * np.linalg.norm(self.r)))

        if np.dot(self.r, self.v) < 0:
            self.nu_0 = 2 * math.pi - self.nu_0

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

    def _find_perigee_or_impact(self):

        # solving for position vector at perigee or impact
        self.impact = False
        self.perigee_or_impact_nu = 0
        self.perigee_radius = ( self.a * ( 1 - np.linalg.norm(self.e)**2 ) ) / ( 1 + np.linalg.norm(self.e) * np.cos(self.perigee_or_impact_nu))

        self.r_p_or_impact_perifocal = np.array([self.perigee_radius*np.cos(self.perigee_or_impact_nu), 
                                                 self.perigee_radius*np.sin(self.perigee_or_impact_nu), 
                                                 0])

        self._solve_for_p_or_impacterifocal_to_ijk_matrix()

        self.r_p_or_impact = np.dot(self.perifocal_to_ijk_matrix, self.r_p_or_impact_perifocal)

        # checking to see if an impact occured
        if np.linalg.norm(self.r_p_or_impact) < 1:
            self.impact = True
            self._solve_for_impact_point()
            return

        # solving for velocity at perigee
        self.v_p_or_impact_perifocal = np.array([math.sqrt(1/self.p) * -1 * np.sin(self.perigee_or_impact_nu),
                                                 np.linalg.norm(self.e) + np.cos(self.perigee_or_impact_nu),
                                                 0])

        self.v_p_or_impact = np.dot(self.perifocal_to_ijk_matrix, self.v_p_or_impact_perifocal) 

    def _solve_for_p_or_impacterifocal_to_ijk_matrix(self):
        
        self.perifocal_to_ijk_matrix = np.array([
                    [
                        np.cos(self.omega) * np.cos(self.w) - np.sin(self.omega) * np.sin(self.w) * np.cos(self.i),
                        -np.cos(self.omega) * np.sin(self.w) - np.sin(self.omega) * np.cos(self.w) * np.cos(self.i),
                        np.sin(self.omega) * np.sin(self.i)
                    ],
                    [
                        np.sin(self.omega) * np.cos(self.w) + np.cos(self.omega) * np.sin(self.w) * np.cos(self.i),
                        -np.sin(self.omega) * np.sin(self.w) + np.cos(self.omega) * np.cos(self.w) * np.cos(self.i),
                        -np.cos(self.omega) * np.sin(self.i)
                    ],
                    [
                        np.sin(self.w) * np.sin(self.i),
                        np.cos(self.w) * np.sin(self.i),
                        np.cos(self.i)
                    ]
                ])

    def _solve_for_impact_point(self):

        self.perigee_or_impact_nu = np.arccos((self.p-1)/np.linalg.norm(self.e))

        self.r_p_or_impact_perifocal = np.array([self.perigee_radius*np.cos(self.perigee_or_impact_nu), 
                                                 self.perigee_radius*np.sin(self.perigee_or_impact_nu), 
                                                 0])
        
        self.r_p_or_impact = np.dot(self.perifocal_to_ijk_matrix, self.r_p_or_impact)

        self.v_p_or_impact_perifocal = np.array([math.sqrt(1/self.p) * -1 * np.sin(self.perigee_or_impact_nu),
                                                 np.linalg.norm(self.e) + np.cos(self.perigee_or_impact_nu),
                                                 0])

        self.v_p_or_impact = np.dot(self.perifocal_to_ijk_matrix, self.v_p_or_impact_perifocal)

    def _find_time_of_flight(self):
        
        if self.trajectory_type == "circular":
            self._time_of_flight_circular()

        elif self.trajectory_type == "elliptical":
            self._time_of_flight_elliptical()

        elif self.trajectory_type == "parabolic":
            self._time_of_flight_parabolic()

        elif self.trajectory_type == "hyperbolic":
            self._time_of_flight_hyperbolic()

        else:
            self.time_of_flight = None

    def _time_of_flight_circular(self):
        pass

    def _time_of_flight_elliptical(self):

        pass

    def _time_of_flight_parabolic(self):

        pass

    def _time_of_flight_hyperbolic(self):

        pass

    def print_results():

        pass


if __name__ == "__main__":

    test_case = 0

    satellite = Satellite( test_case_positions[test_case], test_case_velocities[test_case] )

    print("\ntest_case %d" % (test_case))

    print("satellite.r [ %.4f, %.4f, %.4f ] " % (satellite.r[0], satellite.r[1], satellite.r[2]))
    print("satellite.v [ %.4f, %.4f, %.4f ] " % (satellite.v[0], satellite.v[1], satellite.v[2]))

    print("satellite.r_p_or_impact [ %.4f, %.4f, %.4f ] " % (satellite.r_p_or_impact[0], satellite.r_p_or_impact[1], satellite.r_p_or_impact[2]))
    print("satellite.v_p_or_impact [ %.4f, %.4f, %.4f ] " % (satellite.v_p_or_impact[0], satellite.v_p_or_impact[1], satellite.v_p_or_impact[2]))
    print("satellite trajectory type %s\n\n" % (satellite.trajectory_type))
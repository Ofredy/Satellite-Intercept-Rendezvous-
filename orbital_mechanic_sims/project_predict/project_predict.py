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
import datetime

# Library Imports
import numpy as np

# Our Imports
from project_predict_constants import *


class Satellite:

    def __init__(self, position_vector: list, velocity_vector: list):
        
        self.r_0 = np.array(position_vector)
        self.v_0 = np.array(velocity_vector)

        self._find_orbital_parameters()
        self._find_trajectory_type()
        self._find_perigee_or_impact()
        self._find_time_of_flight()

    def _find_orbital_parameters(self):

        # angular momentum
        self.h = np.cross(self.r_0, self.v_0)

        # semi-latus rectus 
        self.p = np.linalg.norm(self.h)**2

        # eccentricity
        self.e = (np.linalg.norm(self.v_0)**2 - (1/np.linalg.norm(self.r_0))) * self.r_0 - np.dot(self.r_0, self.v_0) * self.v_0

        # semi-major axis -> 1 / a = lamda
        self.lamda = ( 1 - np.linalg.norm(self.e)**2 ) / self.p

        # inclination
        self.i = np.arccos( self.h[k_index] / np.linalg.norm(self.h) )

        # node vector
        self.n = np.cross(k_unit, self.h)

        # longitude of ascending node 
        if abs(np.linalg.norm(self.n)) < 1e-5:
            # ascending node is undefined
            self.omega = 0

        else:
            self.omega = np.arccos( self.n[i_index] / np.linalg.norm(self.n) )

            if self.n[j_index] < 0:
                self.omega = 2 * np.pi - self.omega

        # argument of perigee
        if abs(np.linalg.norm(self.n)) < 1e-5 or abs(np.linalg.norm(self.e)) < 1e-5:
            # argument of perigee undefined
            self.w = 0

        else:
            self.w = np.arccos( np.dot(self.n, self.e) / ( np.linalg.norm(self.n)*np.linalg.norm(self.e) ) )    

            if self.e[k_index] < 0:
                self.w = 2 * np.pi - self.w

        # True anomaly
        if abs(np.linalg.norm(self.e)) < 1e-5:
            # circular orbit
            self.nu_0 = np.dot(self.r_0, self.v_0) / (np.linalg.norm(self.r_0)*np.linalg.norm(self.v_0))

        else:
            self.nu_0 = np.arccos( np.dot(self.e, self.r_0) / ( np.linalg.norm(self.e)*np.linalg.norm(self.r_0) ) )

        if np.dot(self.r_0, self.v_0) < 0:
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
        self.nu = 0
        self.perigee_radius = ( 1/self.lamda ) * ( 1 - np.linalg.norm(self.e) )

        self.r_p_perifocal = np.array([self.perigee_radius*np.cos(self.nu), 
                                       0, 
                                       0])

        self._solve_for_perifocal_to_ijk_matrix()

        self.r = np.dot(self.perifocal_to_ijk_matrix, self.r_p_perifocal)

        # checking to see if an impact occured
        if np.linalg.norm(self.r) < 1:
            self.impact = True
            self._solve_for_impact_point()
            return

        # solving for velocity at perigee
        self.v_p_perifocal = np.sqrt( 1/self.p ) * np.array([0,
                                                             np.linalg.norm(self.e) + np.cos(self.nu),
                                                             0])

        self.v = np.dot(self.perifocal_to_ijk_matrix, self.v_p_perifocal) 

    def _solve_for_perifocal_to_ijk_matrix(self):
        
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

        self.nu = np.arccos( (self.p-1) / np.linalg.norm(self.e) )

        self.r_impact_perifocal = np.array([np.cos(self.nu), 
                                            np.sin(self.nu), 
                                            0])
        
        self.r = np.dot(self.perifocal_to_ijk_matrix, self.r_impact_perifocal)

        self.v_impact_perifocal = math.sqrt(1/self.p) * np.array([-1 * np.sin(self.nu),
                                                                  np.linalg.norm(self.e) + np.cos(self.nu),
                                                                  0])

        self.v = np.dot(self.perifocal_to_ijk_matrix, self.v_impact_perifocal)

    def _find_time_of_flight(self):
        
        if self.nu_0 > self.nu:
            self.d_nu = 2*math.pi + abs(self.nu - self.nu_0)

        else:
            self.d_nu = abs(self.nu - self.nu_0)

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
        
        if self.nu_0 > self.nu:
            self.time_of_flight = ( 2*math.pi - ( self.nu_0 - self.nu ) ) / ( math.sqrt( 1 / (1/self.lamda)**3 ) )

        else:
            self.time_of_flight = ( self.nu - self.nu_0 ) / ( math.sqrt( 1 / (1/self.lamda)**3 ) )

    def _time_of_flight_elliptical(self):

        # calculating eccentric anomalies & change in true anomaly
        self.E_0 = np.arccos( ( np.linalg.norm(self.e) + np.cos(self.nu_0) ) / ( 1 + np.linalg.norm(self.e) * np.cos(self.nu_0) ) )
        self.E = np.arccos( ( np.linalg.norm(self.e) + np.cos(self.nu) ) / ( 1 + np.linalg.norm(self.e) * np.cos(self.nu) ) )

        if self.nu_0 > self.nu:
            self.time_of_flight = math.sqrt( (1/self.lamda)**3 ) * ( 2*math.pi + ( self.E - np.linalg.norm(self.e)*np.sin(self.E) ) - ( self.E_0 - np.linalg.norm(self.e)*np.sin(self.E_0) ) )  

        else:
            self.time_of_flight = math.sqrt( (1/self.lamda)**3 ) * ( ( self.E - np.linalg.norm(self.e)*np.sin(self.E) ) - ( self.E_0 - np.linalg.norm(self.e)*np.sin(self.E_0) ) )

    def _time_of_flight_parabolic(self):

        # calculating parabolic eccentric anomalies & change in true anomaly        
        self.D_0 = math.sqrt( self.p ) * np.tan( self.nu_0 / 2 )
        self.D = math.sqrt( self.p ) * np.tan( self.nu / 2 )

        self.time_of_flight = ( 1 / 2 ) * ( ( self.p*self.D + (1/3)*self.D**3 ) - ( self.p*self.D_0 + (1/3)*self.D_0**3 ) )

    def _time_of_flight_hyperbolic(self):

        # calculating hyperbolic eccentric anomalies & change in true anomaly 
        self.F_0 = np.arccosh( ( np.linalg.norm(self.e) + np.cos(self.nu_0) ) / ( 1 + np.linalg.norm(self.e)*np.cos(self.nu_0) ) )
        self.F = np.arccosh( ( np.linalg.norm(self.e) + np.cos(self.nu) ) / ( 1 + np.linalg.norm(self.e)*np.cos(self.nu) ) )

        self.time_of_flight = math.sqrt( ( -1 * (1/self.lamda) )**3 ) * ( ( np.linalg.norm(self.e)*np.sinh(self.F)-self.F ) - ( np.linalg.norm(self.e)*np.sinh(self.F_0)-self.F_0 ) )

    def _energy_check(self):

        # only used for debugging

        # total specific orbital energy comparison
        self.SE_0 = (np.linalg.norm(self.v_0)**2 / 2) - (1 / np.linalg.norm(self.r_0)) 
        self.SE = (np.linalg.norm(self.v)**2 / 2) - ( 1 / np.linalg.norm(self.r))

        print("E_0 = %f, E = %f" % (self.SE_0, self.SE))

        # angular momentum comparison
        self.h_0 = np.cross(self.r_0, self.v_0)
        self.h = np.cross(self.r, self.v)

        print("h_0 = %f, h = %f" % ( np.linalg.norm(self.h_0), np.linalg.norm(self.h) ))

    def print_results(self):

        print("\ntest_case %d" % (test_case))
        print("satellite.r [ %.4f, %.4f, %.4f ] " % (self.r_0[0], self.r_0[1], self.r_0[2]))
        print("satellite.v [ %.4f, %.4f, %.4f ] " % (self.v_0[0], self.v_0[1], self.v_0[2]))

        print("\ntrajectory type %s" % (self.trajectory_type))

        if self.impact:
            print("impact occured")
            
            print("position at impact [ %.4f, %.4f, %.4f ]" % (self.r[0], self.r[1], self.r[2]))
            print("velocity at impact [ %.4f, %.4f, %.4f ]" % (self.v[0], self.v[1], self.v[2]))
            
            print("time till impact is %s" % (datetime.timedelta(seconds=abs(self.time_of_flight)*TU_TO_SEC)))
            print("change in true anomaly is %.4f" % (math.degrees(self.d_nu)))
        
        else:
            print("position at perigee [ %.4f, %.4f, %.4f ]" % (self.r[0], self.r[1], self.r[2]))
            print("velocity at perigee [ %.4f, %.4f, %.4f ]" % (self.v[0], self.v[1], self.v[2]))

            print("time till perigee is %s" % (datetime.timedelta(seconds=self.time_of_flight*TU_TO_SEC)))
            print("change in true anomaly is %.4f" % (math.degrees(self.d_nu)))


if __name__ == "__main__":

    test_case = 0

    satellite = Satellite( test_case_positions[test_case], test_case_velocities[test_case] )
    satellite.print_results()

    satellite._energy_check()

    import pdb; pdb.set_trace()

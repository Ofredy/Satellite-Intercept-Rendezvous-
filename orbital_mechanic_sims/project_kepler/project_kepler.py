####################### File Description #######################
# Given the position and velocity vectors of a satellite at a particular
# instant of time, you are to determine the position and velocity vectors
# after an interval of time, t


# System imports
import sys
import math

# Library imports
import numpy as np

# Our imports
from project_kepler_constants import *


class Kepler():

    def __init__(self, r_0, v_0, t):
        
        self.r_0 = np.array(r_0)
        self.v_0 = np.array(v_0)
        self.t = t

        self._find_orbital_parameters()
        self._find_trajectory_type()

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

    def _init_x(self):

        self.r_0_norm = np.linalg.norm(self.r_0) 

        # self.lamda = 1 / a -> to avoid dividing by 0
        self.lamda = ( 2 / self.r_0_norm - np.linalg.norm(self.v_0)**2 )

        # solving for the initial guess of x
        if self.trajectory_type == "circular":
            self.x_n = self.nu_0

        elif self.trajectory_type == "elliptical":
            self.x_n = self.t * self.lamda

        elif self.trajectory_type == "parabolic":
            self.x_n = np.sqrt(self.r_0_norm)

        elif self.trajectory_type == "hyperbolic":
            self.x_n = np.sign(self.t) * np.sqrt(-1 * (1/self.lamda)) * np.log( (-2*self.t) / ( (1/self.lamda) * ( np.dot(self.r_0, self.v_0) + np.sign(self.t)*np.sqrt(-1*(1/self.lamda))*(1 - self.r_0_norm*self.lamda)  ) ) )

        else:
            sys.exit("trajectory type not supported")

    def _update_c_and_s(self):

        # updating z
        self.z_n = self.x_n**2 * self.lamda

        if abs(self.z_n) < 1e-5:
            # z close to 0
            self.c_n = (1/math.factorial(2)) - (self.z_n/math.factorial(4)) + (self.z_n**2/math.factorial(6))
            self.s_n = (1/math.factorial(3)) - (self.z_n/math.factorial(5)) + (self.z_n**2/math.factorial(7))

        elif self.z_n > 0:
            # positive z_n
            self.c_n = ( 1 - np.cos(np.sqrt(self.z_n)) ) / self.z_n
            self.s_n = ( np.sqrt(self.z_n) - np.sin(np.sqrt(self.z_n)) ) / np.sqrt(self.z_n**3)

        elif self.z_n < 0:
            # negative z_n
            self.c_n = ( 1 - np.cosh(np.sqrt(-1 * self.z_n)) ) / self.z_n
            self.s_n = ( np.sinh(np.sqrt(-1*self.z_n)) - np.sqrt(-1*self.z_n) ) / np.sqrt( (-1*self.z_n)**3 )

        else:
            sys.exit("error in creating c and z")

    def _solve_for_t_n(self):

        self.t_n = np.dot(self.r_0, self.v_0) * self.x_n**2 * self.c_n + ( 1 - self.r_0_norm*self.lamda )*self.x_n**3*self.s_n + self.r_0_norm*self.x_n 

    def _calculate_t_error(self):

        self.t_error = abs( (self.t - self.t_n) / self.t )

    def _solve_for_dt_n(self):

        # dt_n is short for dt_n/dx_n
        self.dt_n = self.x_n**2*self.c_n + np.dot(self.r_0, self.v_0)*self.x_n*( 1 - self.z_n*self.s_n ) + self.r_0_norm*( 1 - self.z_n*self.c_n )

    def _update_x(self):

        self.x_n += ( self.t - self.t_n ) / self.dt_n

    def _solve_for_r_and_v(self):

        # solving for r
        self.f = 1 - ( self.x_n**2/self.r_0_norm )*self.c_n
        self.g = self.t - self.x_n**3*self.s_n

        self.r = self.f * self.r_0 + self.g * self.v_0

        # solving for v
        self.f_p = ( 1 / (self.r_0_norm*np.linalg.norm(self.r)) ) * self.x_n * ( self.z_n*self.s_n - 1)
        self.g_p = 1 - (self.x_n**2/np.linalg.norm(self.r)) * self.c_n

        self.v = self.f_p * self.r_0 + self.g_p * self.v_0

    def _step_summary(self):

        # only used for debugging
        print("x_%d = %.5f, t_%d = %.5f, dt_%d = %.5f" % (self.counter, self.x_n, self.counter, self.t_n, self.counter, self.dt_n))

    def _f_and_g_check(self):

        # only used for debugging
        print("this should be equal to 1: f*g_p - f_p*g = %.5f" % ( self.f*self.g_p - self.f_p*self.g ))

    def _energy_check(self):

        # only used for debugging

        # total specific orbital energy comparison
        self.E_0 = (np.linalg.norm(self.v_0)**2 / 2) - (1 / self.r_0_norm) 
        self.E = (np.linalg.norm(self.v)**2 / 2) - ( 1 / np.linalg.norm(self.r))

        print("E_0 = %f, E = %f" % (self.E_0, self.E))

        # angular momentum comparison
        self.h_0 = np.cross(self.r_0, self.v_0)
        self.h = np.cross(self.r, self.v)

        print("h_0 = %f, h = %f" % ( np.linalg.norm(self.h_0), np.linalg.norm(self.h) ))

    def kepler_problem(self):

        self._init_x()
        self.counter = 0

        while True:

            self._update_c_and_s()
            self._solve_for_t_n()
            self._calculate_t_error()

            if self.t_error <= 10e-7:
                # solution converged
                self._solve_for_r_and_v()

                print("solution converged")
                print("position is [ %.8f, %.8f, %.8f ]" % (self.r[0], self.r[1], self.r[2]))
                print("velocity is [ %.8f, %.8f, %.8f ]" % (self.v[0], self.v[1], self.v[2]))
                self._f_and_g_check()
                self._energy_check()
                return

            elif self.counter >= 50:
                print("kepler problem did not converge")
                return

            self._solve_for_dt_n()
            self._update_x()
            self.counter += 1

            self._step_summary()


if __name__ == "__main__":

    test_case_idx = 1

    satellite = Kepler( test_case_positions[test_case_idx], test_case_velocities[test_case_idx], test_case_t[test_case_idx] )
    satellite.kepler_problem()

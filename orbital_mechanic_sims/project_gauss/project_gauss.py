####################### File Description #######################
# A satellite in a conic orbit leaves position rl and arrives at position 
# r2 at time dt later. It can travel the "short way" (nu < pi') or the "long
# way" (nu >= pi). Given rl r2 , dt, and sign (pi-nu), find v1 & v2 


# System imports
import sys
import math

# Library imports
import numpy as np

# Our imports
from project_gauss_constants import *


class Gauss:

    def __init__(self):
        pass

    def _evaluate_A(self):

        self.d_nu = np.arccos( (np.dot(self.r_1, self.r_2) / (np.linalg.norm(self.r_1)*np.linalg.norm(self.r_2))) )

        if self.dm == -1:
            self.d_nu = 2*np.pi - self.d_nu

        if abs(self.d_nu) < 1e-3:
            self.A = self.dm * np.sqrt( np.linalg.norm(self.r_1)*np.linalg.norm(self.r_2) * ( 1 + np.cos(self.d_nu) ) )

        else:
            self.A = ( np.sqrt(np.linalg.norm(self.r_1)*np.linalg.norm(self.r_2)) * np.sin(self.d_nu) ) / np.sqrt( 1 - np.cos(self.d_nu) )

    def _init_z(self):
        
        self.z_n = np.random.uniform(-0.5, 0.5)

    def _update_c_and_s(self):

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

    def _update_c_p_and_s_p(self):

        if 0 <= self.z_n and abs(self.z_n) < 1e-3:
            self.c_p_n = ( 1 / math.factorial(4) ) + ( (2*self.z_n) / math.factorial(6) ) - (( 3*self.z_n**2 ) / math.factorial(8) ) + ( ( 4*self.z_n**3 ) / math.factorial(10) )
            self.s_p_n = ( 1 / math.factorial(5) ) + ( (2*self.z_n) / math.factorial(7) ) - (( 3*self.z_n**2 ) / math.factorial(9) ) + ( ( 4*self.z_n**3 ) / math.factorial(11) )

        else:
            self.c_p_n = ( 1 / 2*self.z_n) * ( 1 - self.z_n*self.s_n - 2*self.c_n )
            self.s_p_n = ( 1 / 2*self.z_n) * ( self.c_n - 3*self.s_n )

    def _update_y(self):
        
        self.y_n = np.linalg.norm(self.r_1) + np.linalg.norm(self.r_2) - self.A * ( ( 1 - self.z_n*self.s_n ) / np.sqrt(self.c_n) )

    def _update_x(self):
        
        self.x_n = np.sqrt( self.y_n / self.c_n )

    def _solve_for_t_n(self):
        
        self.t_n = self.x_n**3 * self.s_n + self.A * np.sqrt(self.y_n)

    def _solve_for_dt_n(self):
        
        # dt_n is short for dt_n/dz_n
        self.dt_n = self.x_n**3 * ( self.s_p_n - (( 3*self.s_n*self.c_p_n ) / ( 2*self.c_n )) ) +  ( self.A/8 ) * ( (3*self.s_n*np.sqrt(self.y_n))/self.c_n + self.A/self.x_n )

    def _calculate_t_error(self):
        
        if 1 <= self.t and self.t <= 10e6:
            self.t_error = abs( (self.t - self.t_n) / self.t )

        else:
            self.t_error = abs( self.t - self.t_n )

    def _update_z(self):
        
        self.z_n += (self.t - self.t_n) / self.dt_n

    def _solve_for_v1_and_v2(self):

        self.f = 1 - self.y_n / np.linalg.norm(self.r_1)
        self.g = self.A * np.sqrt( self.y_n  )
        self.g_p = 1 - self.y_n / np.linalg.norm(self.r_2)

        self.v_1 = ( self.r_2 - self.f*self.r_1 ) / self.g

        self.v_2 = ( self.g_p*self.r_2 - self.r_1 ) / self.g 

    def _step_summary(self):

        # only used for debugging
        print("z_%d = %.4f, y_%d = %.4f, x_%d = %.4f, t_%d = %.4f, dt_%d = %.4f, t_error_%d = %.4f" % ( self.counter, self.z_n, self.counter, self.y_n, self.counter, self.x_n, self.counter, self.t_n, self.counter, self.dt_n, self.counter, self.t_error ))

    def _energy_check(self):

        # only used for debugging
        
        # total specific orbital energy comparison
        self.E_1 = (np.linalg.norm(self.v_1)**2 / 2) - (1 / np.linalg.norm(self.r_1)) 
        self.E_2 = (np.linalg.norm(self.v_2)**2 / 2) - ( 1 / np.linalg.norm(self.r_2))

        print("E_0 = %f, E = %f" % (self.E_1, self.E_2))

        # angular momentum comparison
        self.h_1 = np.cross(self.r_1, self.v_1)
        self.h_2 = np.cross(self.r_2, self.v_2)

        print("h_0 = %f, h = %f" % ( np.linalg.norm(self.h_1), np.linalg.norm(self.h_2) ))

    def gauss_problem(self, r_1, r_2, t, dm):

        self.r_1 = np.array(r_1)
        self.r_2 = np.array(r_2)
        self.t = t
        self.dm = dm

        self._evaluate_A()
        self._init_z()

        self.counter = 0

        while True:

            self._update_c_and_s()
            self._update_y()
            self._update_x()
            self._solve_for_t_n()
            self._calculate_t_error()

            if self.t_error < 10e-6:
                # solution converged
                self._solve_for_v1_and_v2()

                print("solution converged")
                print("v1 is [ %.4f, %.4f, %.4f ]" % (self.v_1[0], self.v_1[1], self.v_1[2]))
                print("v2 is [ %.4f, %.4f, %.4f ]" % (self.v_2[0], self.v_2[1], self.v_2[2]))

                self._energy_check()

                return

            elif self.counter >= 50:
                print("gauss problem did not converge")
                return

            self._update_c_p_and_s_p()
            self._solve_for_dt_n()

            self._step_summary()

            self._update_z()

            self.counter += 1


if __name__ == "__main__":

    test_case_idx = 0

    gauss = Gauss()
    
    #gauss.gauss_problem(test_case_r1[test_case_idx], test_case_r2[test_case_idx], test_case_t[test_case_idx], test_case_dm[test_case_idx])
    
    gauss.gauss_problem([0.5, 0.6, 0.7], [0, 1, 0], 0.9667663, 1)
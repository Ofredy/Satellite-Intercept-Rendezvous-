####################### File Description #######################
# Given the position and velocity vectors of a satellite at a particular
# instant of time, you are to determine the position and velocity vectors
# after an interval of time, 61.


# System imports
import math

# Library imports
import numpy as np

# Our imports
from project_kepler_constants import *
from ..project_predict.project_predict import Satellite


class Kepler(Satellite):

    def __init__(self, r_0, v_0, dt):
        
        self.r_0 = r_0
        self.v_0 = v_0
        self.dt = dt

        self._find_orbital_parameters()

    
####################### File Description #######################
# The input to TRACK should be
# range, range rate, elevation, elevation rate, azimuth, azimuth rate,
# latitude of the launch site, radius vector of the launch site and local
# sidereal time of the launch site. The output will be the radius and
# velocity vectors of the satellite.


# System imports
import sys
import math

# Library imports
import numpy as np

# Our imports
from project_site_and_track_constants import *


class Track():

    def __init__(self):
        pass
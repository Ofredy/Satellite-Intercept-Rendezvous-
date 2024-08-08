####################### File Description #######################
# Write a computer program that takes as its input radar tracking data
# on a target satellite, location of the tracking site, time of radar observation,
# and the location of an interceptor launch site. The output is to be the
# impulsive velocity change required for both intercept and interceptplus-rendezvous for various combinations of launch time and
# interceptor time-of-flight. Neglect the atmosphere and assume impulsive
# velocity change from the launch site and at the target


# Library imports
import numpy as np

# Our imports
from project_gauss.project_gauss import Gauss
from project_kepler.project_kepler import Kepler
from project_site_and_track.project_site import site
from project_site_and_track.project_track import track
####################### File Description #######################
# The input to SITE should be latitude,
# altitude and local sidereal time of the launch site and the output should
# be position and velocity of the site.


# System imports
import sys
import math

# Library imports
import numpy as np

# Our imports
from project_site_and_track_constants import *


class Site():

    def __init__(self):
        pass

    def _compute_local_side_real_time(self, day, universal_time):

        """
        Computes the Local Sidereal Time (LST).     
        Parameters:
        longitude (float): Longitude of the launch site in radians (positive to the east).
        day (int): Day of the year where 1 January is day zero.
        ut (float): Universal Time in the format HHMM.SS.       
        Returns:
        float: Local Sidereal Time in radians.
        """

        PI = np.pi
        GST = 1.74933340        

        # Parse UT from the format "HHMM:SS"
        hours = int(universal_time[:2])
        minutes = int(universal_time[2:4])
        seconds = float(universal_time[5:])

        # Calculate the fractional day D
        fractional_day = day + hours / 24 + minutes / 1440 + seconds / 864

        # Calculate Local Sidereal Time
        self.local_side_real_time = self.longitude + GST + 1.0027379093 * 2 * PI * fractional_day  

        self.local_side_real_time = self.local_side_real_time % ( 2 * np.pi ) 

    def site_position_and_velocity(self, latitude_deg, longitude_deg, altitude_ft, day, universal_time):

        self.latitude = np.deg2rad(latitude_deg)
        self.longitude = np.deg2rad(longitude_deg)

        self.altitude = altitude_ft * FT_TO_CANONICAL

        self.x_site = ( 1 / ( np.sqrt( 1 - (EARTH_ECCENTRICITY**2) * (np.sin(self.latitude)**2) ) ) )  + self.altitude
        self.z_site = ( 1 - EARTH_ECCENTRICITY**2 ) / (  np.sqrt( 1 - (EARTH_ECCENTRICITY**2) * (np.sin(self.latitude)**2) ) ) + self.altitude

        self._compute_local_side_real_time(day, universal_time)

        self.r_site = np.array([
                                    self.x_site * np.cos(self.latitude) * np.cos(self.local_side_real_time),
                                    self.x_site * np.cos(self.latitude) * np.sin(self.local_side_real_time),
                                    self.z_site * np.sin(self.latitude)
                                ])

        self.v_site = np.cross(EARTH_ANGULAR_VELOCITY, self.r_site)


if __name__ == "__main__":

    site = Site()

    site.site_position_and_velocity(39.007, -104.883, 7180, 244, "0317:02")

    #site.site_position_and_velocity(0, -57.296, 20901.33, 1, "0600:00")

    print("r_site is [ %.4f, %.4f, %.4f ]" % (site.r_site[0], site.r_site[1], site.r_site[2]))
    print("v_site is [ %.4f, %.4f, %.4f ]" % (site.v_site[0], site.v_site[1], site.v_site[2]))    

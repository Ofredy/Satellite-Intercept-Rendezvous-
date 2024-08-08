####################### File Description #######################
# The input to TRACK should be
# range, range rate, elevation, elevation rate, azimuth, azimuth rate,
# latitude of the launch site, radius vector of the launch site and local
# sidereal time of the launch site. The output will be the radius and
# velocity vectors of the satellite.


# Library imports
import numpy as np

# Our imports
from project_site_and_track_constants import *
from project_site import Site


class Track():

    def __init__(self):
        pass

    def _SEZ_to_IJK_matrix(self):

        self.sez_to_ijk_matrix = np.array( [ [ np.sin(self.site_latitude) * np.cos(self.site_local_sidereal_time), -np.sin(self.site_local_sidereal_time), np.cos(self.site_latitude) * np.cos(self.site_local_sidereal_time) ],
                                             [ np.sin(self.site_latitude) * np.sin(self.site_local_sidereal_time), np.cos(self.site_local_sidereal_time), np.cos(self.site_latitude) * np.sin(self.site_local_sidereal_time) ],
                                             [ -np.cos(self.site_latitude), 0, np.sin(self.site_latitude) ] ] )

    def track(self, range, range_p, elevation, elevation_p, azimuth, azimuth_p, r_site, site_latitude, site_local_sidereal_time):

        self.site_latitude = site_latitude
        self.site_local_sidereal_time = site_local_sidereal_time

        range *= KM_TO_CANONICAL
        range_p *= KM_TO_CANONICAL * SEC_TO_TU

        elevation = np.deg2rad(elevation)
        elevation_p = np.deg2rad(elevation_p) * SEC_TO_TU

        azimuth = np.deg2rad(azimuth)
        azimuth_p = np.deg2rad(azimuth_p) * SEC_TO_TU

        self.r_satellite_rel_site_sez = np.array([  -range * np.cos(elevation) * np.cos(azimuth),
                                                    range * np.cos(elevation) * np.sin(azimuth),
                                                    range * np.sin(elevation) ])

        self.v_satellite_rel_site_sez = np.array([  -range_p * np.cos(elevation) * np.cos(azimuth) 
                                                    + range * np.sin(elevation) * elevation_p * np.cos(azimuth) 
                                                    + range * np.cos(elevation) * np.sin(azimuth) * azimuth_p,
                                                    range_p * np.cos(elevation) * np.sin(azimuth)
                                                    - range * np.sin(elevation) * elevation_p * np.sin(azimuth)
                                                    + range * np.cos(elevation) * np.cos(azimuth) * azimuth_p,
                                                    range_p * np.sin(elevation) + range * np.cos(elevation) * elevation_p ])

        self._SEZ_to_IJK_matrix()
        self.r_satellite_rel_site = np.dot(self.sez_to_ijk_matrix, self.r_satellite_rel_site_sez)
        self.v_satellite_rel_site = np.dot(self.sez_to_ijk_matrix, self.v_satellite_rel_site_sez)

        self.r_satellite = r_site + self.r_satellite_rel_site

        self.v_satellite = self.v_satellite_rel_site + np.cross(EARTH_ANGULAR_VELOCITY, self.r_satellite)


if __name__ == "__main__":

    site = Site()

    site.site_position_and_velocity(39.007, -104.883, 7180, 244, "0317:02")

    track = Track()

    r_site = np.array([.20457216, -.75100391 , .62624920])

    track.track(504.68, 2.08, 30.7, 0.07, 105.7, 0.05, r_site, site.latitude, site.local_side_real_time)

    print("r_satellite is [ %.4f, %.4f, %.4f ]" % (track.r_satellite[0], track.r_satellite[1], track.r_satellite[2]))
    print("v_satellite is [ %.4f, %.4f, %.4f ]" % (track.v_satellite[0], track.v_satellite[1], track.v_satellite[2]))    

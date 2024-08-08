####################### File Description #######################
# Write a computer program that takes as its input radar tracking data
# on a target satellite, location of the tracking site, time of radar observation,
# and the location of an interceptor launch site. The output is to be the
# impulsive velocity change required for both intercept and interceptplus-rendezvous for various combinations of launch time and
# interceptor time-of-flight. Neglect the atmosphere and assume impulsive
# velocity change from the launch site and at the target


# Our imports
from project_site_and_track.project_site import Site
from project_site_and_track.project_track import Track
from project_kepler.project_kepler import Kepler
from project_gauss.project_gauss import Gauss
from project_intercept_constants import *


def update_universal_time(initial_universal_time, minutes_passed):

     # Parse the input time string
    hours = int(initial_universal_time[:2])
    minutes = int(initial_universal_time[2:4])
    seconds = int(initial_universal_time[5:7])

    # Calculate total minutes
    total_minutes = hours * 60 + minutes + minutes_passed

    # Calculate the new hours, minutes, and seconds
    new_hours = (total_minutes // 60) % 24  # Modulo 24 to wrap around after 24 hours
    new_minutes = total_minutes % 60

    # Format the updated time string
    new_univerval_time = f"{new_hours:02d}{new_minutes:02d}:{seconds:02d}"

    return new_univerval_time

def calculate_delta_v(v1, v2):

    return np.linalg.norm(v2- v1)

def rendezvous_analysis_summary(launch_time, tof, site_v, satellite_v, direct_v0, direct_v1, retro_v0, retro_v1):

    print("launch_time: %d, tof: %d" % (launch_time, tof))

    if isinstance(direct_v0, np.ndarray) and isinstance(direct_v1, np.ndarray):
        delta_v_site = calculate_delta_v(site_v, direct_v0)
        delta_v_rendezvous = calculate_delta_v(direct_spacecraft_v_1, satellite_v)
        print("direct: delta_v_site: %.5f, delta_v_rendezvous: %.5f" % (delta_v_site, delta_v_rendezvous))

    else:
        print("direct orbit solution did not converge")

    if isinstance(retro_v0, np.ndarray) and isinstance(retro_v1, np.ndarray):
        delta_v_site = calculate_delta_v(site_v, retro_v0)
        delta_v_rendezvous = calculate_delta_v(retro_v1, satellite_v)
        print("retro: delta_v_site: %.5f, delta_v_rendezvous: %.5f" % (delta_v_site, delta_v_rendezvous))

    else:
        print("direct orbit solution did not converge")


if __name__ == "__main__":

    radar_site = Site()
    launch_site = Site()
    track = Track()
    
    kepler = Kepler()
    gauss = Gauss()


    # Solve for satellite position and velocity
    radar_site_r, radar_site_v = radar_site.site_position_and_velocity(RADAR_TRACKING_SITE["latitude"], RADAR_TRACKING_SITE["longitude"], RADAR_TRACKING_SITE["altitude"], 
                                                                 RADAR_TRACKING_SITE["day"], RADAR_TRACKING_SITE["universal_time"])

    satellite_r_0, satellite_v_0 = track.track(RADAR_MEASUREMENTS["range"], RADAR_MEASUREMENTS["range_p"], RADAR_MEASUREMENTS["elevation"], RADAR_MEASUREMENTS["elevation_p"],
                                                RADAR_MEASUREMENTS["azimuth"], RADAR_MEASUREMENTS["azimuth_p"], radar_site_r, radar_site.latitude, radar_site.local_side_real_time)

    # Begining launch time & time of flight analysis
    for launch_time in range(FIRST_REC_T, FIRST_REC_T + REC_T_INCREMENT_SIZE*NUM_REC_TIMES + REC_T_INCREMENT_SIZE, REC_T_INCREMENT_SIZE):

        for tof in range(FIRST_TOF, FIRST_TOF + TOF_INCREMENT_SIZE*NUM_TOFS + TOF_INCREMENT_SIZE, TOF_INCREMENT_SIZE):

            # solving for launch site position & velocity at launch time
            launch_universal_time = update_universal_time(RADAR_TRACKING_SITE["universal_time"], launch_time)
            launch_site_r, launch_site_v = launch_site.site_position_and_velocity(LAUNCH_SITE["latitude"], LAUNCH_SITE["longitude"], LAUNCH_SITE["altitude"], RADAR_TRACKING_SITE["day"], launch_universal_time)

            # solving for satellite position & velocity after time = launch_time + tof
            t = launch_time + tof
            satellite_r_1, satellite_v_1 = kepler.kepler_problem(satellite_r_0, satellite_v_0, t)

            # solving for the needed velocities of our spacecraft

            # direct orbit
            direct_spacecraft_v_0, direct_spacecraft_v_1 = gauss.gauss_problem(launch_site_r, satellite_r_1, t, 1)

            # retrograde orbit
            retro__spacecraft_v_0, retro_spacecraft_v_1 = gauss.gauss_problem(launch_site_r, satellite_r_1, t, -1)

            # finding delta_vs
            rendezvous_analysis_summary(launch_time, tof, launch_site_v, satellite_v_1, direct_spacecraft_v_0, direct_spacecraft_v_1, retro__spacecraft_v_0, retro_spacecraft_v_1)

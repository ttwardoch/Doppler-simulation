import numpy as np
import matplotlib.pyplot as plt

G = 6.6743 * 10**(-11)
c = 300000000
EARTH_MASS = 5.97219 * 10**24
EARTH_RADIUS = 6371000  # in m
EARTH_SPEED = 0#2*np.pi*EARTH_RADIUS / 24 / 3600  # speed of the equator in m/s
SATELLITE_HEIGHT = 600000  # in m
SATELLITE_RADIUS = EARTH_RADIUS + SATELLITE_HEIGHT
ORBITAL_SPEED = np.sqrt(G*EARTH_MASS/SATELLITE_RADIUS)  # m/s
ANGULAR_SPEED = np.sqrt(G*EARTH_MASS/SATELLITE_RADIUS) / SATELLITE_RADIUS
FREQUENCY = 1000000000

print(f"Orbital speed: {ORBITAL_SPEED}")
print(f"Earth speed: {EARTH_SPEED}")
print(np.sqrt(EARTH_SPEED**2+ORBITAL_SPEED**2))
print(40046*ANGULAR_SPEED*180/np.pi)

LATITUDE = 0

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)
def cos_of_angle_between(v1, v2):
    """ Returns the cos of angle between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)

def doppler_shift(satellite_angle_perp, satellite_angle_tang):
    """
    :param satellite_angle_perp: inclination of the satellite in degrees
    :param satellite_angle_tang: angular distance of the satellite from the highest point above the horizon (- is before achieving the highest point) in degrees
    :return: value of Doppler shift in Hz
    """
    position_of_satellite = np.array([0, SATELLITE_RADIUS*np.cos(satellite_angle_tang*np.pi/180), SATELLITE_RADIUS*np.sin(satellite_angle_tang*np.pi/180)])
    velocity_of_satellite = np.array([0, ORBITAL_SPEED*np.sin(-satellite_angle_tang*np.pi/180), ORBITAL_SPEED*np.cos(satellite_angle_tang*np.pi/180)])
    position_of_UE = np.array([EARTH_RADIUS * np.sin(satellite_angle_perp*np.pi/180), EARTH_RADIUS*np.cos(satellite_angle_perp*np.pi/180), 0])
    velocity_of_UE = np.array([EARTH_SPEED*np.cos(satellite_angle_perp*np.pi/180), EARTH_SPEED * np.sin(-satellite_angle_perp*np.pi/180),0])
    position_of_UE_minus_satellite = position_of_UE - position_of_satellite
    velocity_of_satellite_minus_UE = velocity_of_satellite - velocity_of_UE

    # Check if above horizon
    if np.dot((position_of_satellite - position_of_UE), position_of_UE) < 0:
        raise ValueError

    doppler = FREQUENCY / c * np.dot(velocity_of_satellite_minus_UE, unit_vector(position_of_UE_minus_satellite))
    return doppler

for inclination in [0, 2, 4, 6, 8, 10, 15, 20]:
    ts = list(np.linspace(-3000, 3000, num=10000))
    dopplers = []
    tsbis = []
    for t in ts:
        try:
            dopplers.append(doppler_shift(inclination, t*ANGULAR_SPEED*360/2/np.pi))
            tsbis.append(t)
        except ValueError:
            pass
    plt.plot(tsbis, dopplers, label=f"{inclination}°")
plt.xlabel("Time[s]")
plt.ylabel("Doppler shift per GHz [kHz/Ghz]")
plt.title("Doppler shift vs time for a satellite above the horizon")
plt.legend()
plt.tight_layout()
plt.show()
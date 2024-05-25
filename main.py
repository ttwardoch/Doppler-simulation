import numpy as np
import matplotlib.pyplot as plt

EARTH_RADIUS = 6371000 # in m
SATELLITE_HEIGHT = 600000 # in m
SATELLITE_RADIUS = EARTH_RADIUS + SATELLITE_HEIGHT

G = 6.6743 * 10**(-11)
EARTH_MASS = 5.97219 * 10**24
FREQUENCY = 1 # GHz
dt = 1 # seconds
iterations = 10000

ORBITAL_SPEED = np.sqrt(G*EARTH_MASS/SATELLITE_RADIUS) # m/s
print(f"Orbital speed: {ORBITAL_SPEED/1000} km/s")
ANGULAR_SPEED = np.sqrt(G*EARTH_MASS/SATELLITE_RADIUS) / SATELLITE_RADIUS


FREQUENCY *= 1000000000  # convert frequency
ts = np.array(range(iterations)) * dt

theta = 0
for theta in [0, 5, 10, 15, 20]:
    UE_coord = [0, EARTH_RADIUS*np.cos(theta/180 * np.pi), EARTH_RADIUS*np.sin(theta/180 * np.pi)]

    satellite_coords = []
    for t in ts:
        satellite_coords.append([np.sin(t*ANGULAR_SPEED)*SATELLITE_RADIUS, np.cos(t*ANGULAR_SPEED)*SATELLITE_RADIUS, 0])

    satellite_coords = np.array(satellite_coords)
    distance = [np.sqrt(coord[0]**2+(coord[1]-EARTH_RADIUS*np.cos(theta/180 * np.pi))**2+(EARTH_RADIUS*np.sin(theta/180 * np.pi))**2) if np.dot(UE_coord, [coord[0], coord[1], 0]) >= EARTH_RADIUS**2 else 0 for coord in satellite_coords]
    relative_velocity = []
    for i in range(iterations-1):
        speed = (distance[i+1]-distance[i])/dt
        if speed < 10000 and speed > -100000:
            relative_velocity.append(speed)
        else:
            relative_velocity.append(0)

    tss = []
    relative_velocitys = []
    for t, vel in zip(ts[:-1]-6297+505, relative_velocity):
        if vel != 0 and t>-3000:
            tss.append(t)
            relative_velocitys.append(vel)
    tss = np.array(tss)
    plt.plot(tss, -np.array(relative_velocitys)/300000000*FREQUENCY/1000, label=f"{theta}°")
plt.xlabel("Time [s]")
plt.ylabel("Doppler shift per GHz [kHz/GHz]")
plt.legend()
plt.title("Doppler shift for a satellite above the horizon")
plt.tight_layout()
plt.show()


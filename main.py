"""
Main file which executes the code. Runs a simulation for N asteroids under the
gravitational effects of the Sun and Jupiter.
"""
import kirkwoods
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import random
from scipy.integrate import odeint


# !Want to work in units of AU/SolarMass/Years!
## GM = 4pi^2
GM = 4*np.pi**2
MSOL = 3.33e5

test = kirkwoods.Simulation(250, 120, 0.002)
test.run_N_body_sim(test.Sun, test.Jupiter, test.asteroids)

# 2D plotting stuff below
#plt.plot(*test.Jupiter.pos[:2], label="Jupiter")
#plt.plot(*test.Sun.pos[:2], lw=4, label="Sun")
#for asteroid in test.asteroids:
#    plt.plot(*asteroid.pos[:2], ls="dashed")  #, label="Asteroid"
#plt.legend(fontsize=12, frameon=True, fancybox=True, edgecolor="#00AA00", loc="lower right")
#plt.axis([-5.5, 5.5, -5.5, 5.5])
#plt.show()

print len(test.asteroids)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(*test.Sun.pos, lw=4, label="Sun")
for asteroid in test.asteroids:
    ax.plot(*asteroid.pos, ls="dashed")  #, label="Asteroid"
ax.plot(*test.Jupiter.pos, label="Jupiter", c="black", lw=2)
ax.legend(fontsize=12, frameon=True, fancybox=True, edgecolor="#00AA00", loc="lower right")
ax.set_xlim3d(-5.5, 5.5)
ax.set_xlabel("X (AU)")
ax.set_ylim3d(-5.5, 5.5)
ax.set_ylabel("Y (AU)")
ax.set_zlim3d(-.1, .1)  # Unsure what the best values are here, this seems pretty ok
ax.set_zlabel("Z (AU)")

plt.show()

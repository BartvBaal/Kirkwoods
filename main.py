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

test = kirkwoods.Simulation(500, 100, 0.002)
test.run_N_body_sim(test.Sun, test.Jupiter, test.asteroids)
test.visualize()

# 2D plotting stuff below
#plt.plot(*test.Jupiter.pos[:2], label="Jupiter")
#plt.plot(*test.Sun.pos[:2], lw=4, label="Sun")
#for asteroid in test.asteroids:
#    plt.plot(*asteroid.pos[:2], ls="dashed")  #, label="Asteroid"
#plt.legend(fontsize=12, frameon=True, fancybox=True, edgecolor="#00AA00", loc="lower right")
#plt.axis([-5.5, 5.5, -5.5, 5.5])
#plt.show()



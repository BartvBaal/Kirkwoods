"""
Main file which executes the code. Currently does a 2D jupiter-sun simulation
"""
import kirkwoods
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import random
from scipy.integrate import odeint

# !Want to work in units of AU/SolarMass/Years!
## GM = 4pi^2
GM = 4*np.pi**2
MSOL = 3.33e5

test = kirkwoods.Simulation(4, 12, 0.001)
test.run_N_body_sim(test.Sun, test.Jupiter, test.asteroids)

plt.plot(*test.Jupiter.pos, label="Jupiter")
for asteroid in test.asteroids:
    plt.plot(*asteroid.pos, ls="dashed")  #, label="Asteroid"
plt.plot(*test.Sun.pos, lw=4, label="Sun")
plt.legend(fontsize=12, frameon=True, fancybox=True, edgecolor="#00AA00", loc="lower right")
plt.axis([-5.5, 5.5, -5.5, 5.5])
plt.show()

print len(test.Sun.pos[0]), len(test.Jupiter.pos[0]), len(test.asteroids[0].pos[0])

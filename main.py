"""
Main file which executes the code. Currently does a 2D jupiter-sun simulation
"""
import kirkwoods
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import random
from scipy.integrate import odeint

# Probably want to do this initialization as a whole group
# So we can set the total time to be equal for all objects
# !Want to work in units of AU/SolarMass/Years!
## GM = 4pi^2
GM = 4*np.pi**2
MSOL = 3.33e5

test = kirkwoods.Simulation(2, 12, 0.001)
test.run_two_body_sim(test.Sun, test.Jupiter)
for asteroid in test.asteroids:
    test.run_two_body_sim(test.Sun, asteroid)
    plt.plot(*asteroid.pos, label="Asteroid")

#print len(test.Sun.pos[0]), test.Sun.time_step
plt.plot(*test.Jupiter.pos, label="Jupiter")
plt.plot(*test.Sun.pos, lw=8)  # Sun cuz it doesnt show stuff for test.Sun.pos ??
plt.legend(fontsize=12, frameon=True, fancybox=True, edgecolor="#00AA00", loc="lower right")
plt.show()

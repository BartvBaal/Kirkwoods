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
## GM = 4pi^2 => G = 4pi^2
G = 4*np.pi**2
MSOL = 3.33e5

test = kirkwoods.Simulation(0, 15, 0.001)
test.run_two_body_sim(test.Sun, test.Jupiter)

print len(test.Sun.pos[0]), test.Sun.time_step
plt.plot(*test.Jupiter.pos)
plt.plot([0, 0], 'o', lw=8)  # prints something at 1,0 too?????
plt.show()

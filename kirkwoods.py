"""
    Computational Astrophysics: Numerical Integration final assignment
    Contributors:               Bart van Baal, Boris Wolvers
    Filename:                   kirkwoods.py
    Description:
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import random
from scipy.integrate import odeint

# Change border size of plots
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5

class Kirkwoods(object):
    """
    This class contains default initialisations bla bla
    """

    def __init__(self, position, velocity, mass,
                number_of_seconds = 20, time_step = 0.1):
        """
        Initialisations, if none provided: default parameters are used
        """

        # Initial positions
        self.x0 = position[0]
        self.y0 = position[1]

        # Initial velocity
        self.vx0 = velocity[0]
        self.vy0 = velocity[1]

        # Initial mass
        self.mass = mass

        # Duration of the system
        self.number_of_seconds = number_of_seconds

        # The interval (step size or h)
        self.time_step = time_step

        # An array of time intervals (will be same for all of the methods)
        self.time_array = np.arange(0, number_of_seconds, time_step)

    def harm_oscil_analytic(self):
        """
        Solves analytically the harmonic oscillator
        """

        position_analytic = self.amplitude*np.cos(self.time_array)

        return position_analytic

    def harm_oscil_euler_cromer(self):
        """
        Solves ODE using the Euler-Cromer method
        """

        position = [self.x0]
        velocity = [self.v0]

        # Calculates new pos and velo for every time step
        for i in range(len(self.time_array)-1):

            # Calculates new pos and velo using Euler-Cromer method
            vi_plus_1 = velocity[-1] - self.time_step*position[-1]
            xi_plus_1 = position[-1] + self.time_step*vi_plus_1

            position.append(xi_plus_1)
            velocity.append(vi_plus_1)

        return np.asarray(position), np.asarray(velocity)

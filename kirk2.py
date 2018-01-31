"""
    Computational Astrophysics: Numerical Integration final assignment
    Contributors:               Bart van Baal, Boris Wolvers
    Filename:                   kirkwoods.py
    Description:
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import animation
import random
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# Change border size of plots
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5

class Constants(object):
    """
    Defining a class that contains the constants within the problem
    """
    def __init__(self):
        # G is PER SOLAR MASS
        self.GM = 4*np.pi**2
        self.MSOL = 1

        # Jupiter info
        self.smaxis_jup = 5.2044
        self.mass_jup = 1/1047.
        self.ecc_jup = 0.0489

        # Point 0,0 to be the center of mass
        self.offset_cm = self.smaxis_jup / (1+(1/self.mass_jup))
        self.start_vel_jup = math.sqrt(
                ((self.GM)/self.smaxis_jup)*((1+self.ecc_jup)/1-self.ecc_jup))

        # Initial conditions of sun and jupiter (x0, y0, z0) and (vx0, vy0, vz0)
        self.initial_pos_jup = np.asarray([0, self.smaxis_jup*(1-self.ecc_jup)-self.offset_cm, 0])
        self.initial_vel_jup = np.asarray([self.start_vel_jup, 0, 0])

        self.initial_pos_sun = np.asarray([0, -self.offset_cm, 0])
        self.initial_vel_sun = np.asarray([-self.start_vel_jup*self.offset_cm/self.smaxis_jup,0, 0])


class Astro_body(object):
    """
    This class contains default initialisations
    """

    def __init__(self, position, velocity):
        """
        Initialisations, if none provided: default parameters are used
        """

        # Current location and velocity of the objects
        self.pos = np.asarray(position)
        self.vel = np.asarray(velocity)

        # Variable to track the distance to the sun and jupiter
        self.dis_sun = 0
        self.dis_jup = 0

class Kirkwood_solver(object):
    """
    """
    def __init__(self, total_time, time_step, amount_of_asteroids, constants):
        """
        """

        # The interval (step size or h)
        self.time_step = time_step

        # number of iterations
        self.n_iterations = total_time / self.time_step

        self.const = constants

        self.amount_asteroids = amount_of_asteroids

        self.sun_pos = constants.initial_pos_sun
        self.sun_vel = constants.initial_vel_sun
        self.jup_pos = constants.initial_pos_jup
        self.jup_vel = constants.initial_vel_jup

        ast_pos, ast_vel = self.create_asteroids()
        self.asteroids_pos = ast_pos
        self.asteroids_vel = ast_vel

    def create_asteroids(self):
        ast_sema = (2, 5.)

        asteroids_pos = []  # List off asteroids (for now)
        asteroids_vel = []

        # Create all the asteroid Kirkwoods objects
        for amount in range(self.amount_asteroids):
            startecc = 0  # No eccentricity for now
            startloc = random.uniform(*ast_sema)*(1-startecc)
            startvel = np.sqrt(((1+startecc)/(1-startecc))*self.const.GM/startloc)

            # Randomize the starting point of the asteroid's orbit
            orb_loc = random.uniform(0, 2*np.pi)
            startx = np.cos(orb_loc)*startloc
            starty = np.sin(orb_loc)*startloc
            startz = 0  # Start in the jupiter-sun plane
            startvelx = np.sin(orb_loc)*startvel
            startvely = -np.cos(orb_loc)*startvel
            startvelz = random.uniform(-startvel, startvel)*0.01  # Temp range

            asteroids_pos.append([startx, starty, startz])
            asteroids_vel.append([startvelx, startvely, startvelz])

        return np.asarray(asteroids_pos), np.asarray(asteroids_vel)

    def update_planet(self):
        """
        Updates the location and velocity for a planet around the sun for a
        single timestep.
        """
        # Sun and Jupiter should never get a z-location != 0 so only 2D
        distance = math.sqrt(np.sum((self.sun_pos - self.jup_pos)**2))

        self.sun_vel = self.sun_vel - (self.const.GM*self.const.mass_jup*(self.sun_pos - self.jup_pos) / (distance**3))*self.time_step
        self.sun_pos = self.sun_pos + self.sun_vel*self.time_step

        self.jup_vel = self.jup_vel - (self.const.GM*self.const.MSOL*(self.jup_pos - self.sun_pos) / (distance**3))*self.time_step
        self.jup_pos = self.jup_pos + self.jup_vel*self.time_step

    def update_asteroid(self):
        """
        Updates the position for the ast roid *before* the sun and planet have
        been updated. Body1 as star, body2 as planet, body3 as asteroid.
        dimensions should be the same for all objects. Currently works for 3D.
        """
        dis_sun = np.sqrt(np.sum((self.asteroids_pos - self.sun_pos)**2,  axis=1))[:,None]
        dis_jup = np.sqrt(np.sum((self.asteroids_pos - self.jup_pos)**2,  axis=1))[:,None]

        self.asteroids_vel = self.asteroids_vel - self.const.GM*((self.const.MSOL*(self.asteroids_pos - self.sun_pos) / (dis_sun**3)) +
                               (self.const.mass_jup*(self.asteroids_pos - self.jup_pos) / (dis_jup**3)))*self.time_step
        self.asteroids_pos = self.asteroids_pos + self.asteroids_vel*self.time_step

        # Remove the asteroid from the list if it goes too far out
        # Current limit is 6.5AU
        # if body3.dis_sun > 6.5:
        #     self.asteroids.remove(body3)

    def run_N_body_sim(self):
        """
        """

        for i in range(int(self.n_iterations) - 1):
            self.update_asteroid()
            self.update_planet()


c = Constants()
#sun = Astro_body(c.initial_pos_sun, c.initial_vel_sun, c.MSOL, 0)
#jupiter = Astro_body(c.initial_pos_jup, c.initial_vel_jup, c.mass_jup, c.ecc_jup)

#amount_of_asteroids, total_time, time_step)
test = Kirkwood_solver(10, 0.002, 2, c)
test.run_N_body_sim()
print type(test)
#
print "n",test.asteroids_pos

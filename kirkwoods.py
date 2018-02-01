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
        # Units in AU, solar masses and years
        self.GM = 4*np.pi**2
        self.MSOL = 1

        # Jupiter info
        self.smaxis_jup = 5.2044
        self.orbital_period_jup = self.smaxis_jup**(3/2.)
        self.mass_jup = 1/1047.
        self.ecc_jup = 0.0489

        # Point 0,0 to be the center of mass of sun/jupiter system
        self.offset_cm = self.smaxis_jup / (1+(1/self.mass_jup))

        # Initial velocity of jupiter
        self.start_vel_jup = math.sqrt(
                ((self.GM)/self.smaxis_jup)*((1+self.ecc_jup)/1-self.ecc_jup))

        # Initial conditions of jupiter and sun (x0, y0, z0) and (vx0, vy0, vz0)
        self.initial_pos_jup = np.asarray([0,
                            self.smaxis_jup*(1-self.ecc_jup)-self.offset_cm, 0])
        self.initial_vel_jup = np.asarray([self.start_vel_jup, 0, 0])

        self.initial_pos_sun = np.asarray([0, -self.offset_cm, 0])
        self.initial_vel_sun = np.asarray(
                    [-self.start_vel_jup*self.offset_cm/self.smaxis_jup,0, 0])


class Kirkwood_solver(object):
    """
    This class contains important initialisations, creating all of the asteroids
    and solves the orbital motions of jupiter and the asteroids using
    Euler-Cromer
    """
    def __init__(self, total_time, time_step, number_of_asteroids, constants):
        """
        Initialisations
        """

        # The interval (step size or h)
        self.time_step = time_step

        # Number of iterations
        self.n_iterations = total_time / self.time_step

        # Setting the constants of the problem
        self.const = constants

        # Number of asteroids
        self.number_of_asteroids = number_of_asteroids

        # Position and velocities of sun/jupiter (x,y,z), (vx,vy,vz) as numpy
        # arrays, as defined in the constants
        self.sun_pos = constants.initial_pos_sun
        self.sun_vel = constants.initial_vel_sun
        self.jup_pos = constants.initial_pos_jup
        self.jup_vel = constants.initial_vel_jup

        # Obtain all of the asteroids, their positions and velocities
        ast_pos, ast_vel = self.create_asteroids()
        self.asteroids_pos = ast_pos
        self.asteroids_vel = ast_vel


    def create_asteroids(self):
        """
        Creating the asteroids with random intial conditions
        """

        # The range of semi major axis where asteroids should be in
        smaxis_asteroid = (2, 5.)

        # Storing as (x,y,z) and (vx,vy,vz) respectively
        asteroids_pos = []
        asteroids_vel = []

        # Create all the asteroids with varying initial settings
        for amount in range(self.number_of_asteroids):
            startecc = 0  # No eccentricity for now
            startloc = random.uniform(*smaxis_asteroid)*(1-startecc)
            startvel = math.sqrt(((1+startecc)/(1-startecc))*self.const.GM/startloc)

            # Randomize the starting point of the asteroid's orbit
            orb_loc = random.uniform(0, 2*np.pi)
            startx = math.cos(orb_loc)*startloc
            starty = math.sin(orb_loc)*startloc
            startz = 0  # Start in the jupiter-sun plane
            startvelx = math.sin(orb_loc)*startvel
            startvely = -math.cos(orb_loc)*startvel
            startvelz = random.uniform(-startvel, startvel)*0.01  # Temp range

            asteroids_pos.append([startx, starty, startz])
            asteroids_vel.append([startvelx, startvely, startvelz])

        return np.asarray(asteroids_pos), np.asarray(asteroids_vel)


    def update_planet(self):
        """
        Updates the location and velocity for a planet around the sun for a
        single timestep.
        """
        # Determine distance between sun and jupiter
        distance = math.sqrt(np.sum((self.sun_pos - self.jup_pos)**2))

        # Update solar velocity and location; vdot = GM*vector/(distance**3)
        self.sun_vel = self.sun_vel - (self.const.GM*self.const.mass_jup*(self.sun_pos - self.jup_pos) / (distance**3))*self.time_step
        self.sun_pos = self.sun_pos + self.sun_vel*self.time_step

        # Update jupiter velocity and location; vdot = GM*vector/(distance**3)
        self.jup_vel = self.jup_vel - (self.const.GM*self.const.MSOL*(self.jup_pos - self.sun_pos) / (distance**3))*self.time_step
        self.jup_pos = self.jup_pos + self.jup_vel*self.time_step


    def update_asteroid(self):
        """
        Updates the velocities and locations for all asteroids. Should be called
        before the update_planet() function, to prevent off-by-one calculations.
        Will also throw away asteroids which go too far from the sun.
        """
        # Matrix setup for distances to the sun/jupiter for each asteroid
        dis_sun = np.sqrt(np.sum((self.asteroids_pos - self.sun_pos)**2,  axis=1))[:,None]
        dis_jup = np.sqrt(np.sum((self.asteroids_pos - self.jup_pos)**2,  axis=1))[:,None]

        # Check and remove for runaway asteroids ## NOTE: NOT 100% EFFECTIVE
        index_far_asteroids = np.where(dis_sun > 7)[0]
        if len(index_far_asteroids) > 1:
#            print dis_sun[index_far_asteroids], type(index_far_asteroids)
            dis_sun = np.delete(dis_sun, index_far_asteroids, axis=0)
            dis_jup = np.delete(dis_jup, index_far_asteroids, axis=0)

            self.asteroids_pos = np.delete(self.asteroids_pos, index_far_asteroids, axis=0)
            self.asteroids_vel = np.delete(self.asteroids_vel, index_far_asteroids, axis=0)

        # Update asteroids velocities and locations; vdot = GM*vector/(distance**3)
        self.asteroids_vel = (self.asteroids_vel - 
           self.const.GM*((self.const.MSOL*(self.asteroids_pos - self.sun_pos) / (dis_sun**3)) +
           (self.const.mass_jup*(self.asteroids_pos - self.jup_pos) / (dis_jup**3)))*self.time_step)
        self.asteroids_pos = self.asteroids_pos + self.asteroids_vel*self.time_step


    def run_N_body_sim(self):
        """
        Calls the update_asteroid and update_planet functions in order to run
        the simulation for the amount of asteroids specified in the init.
        See those specific functions for more comments on their workings.
        """
        # Perform n_iterations-1 steps (intialization counts for the first step)
        for i in range(int(self.n_iterations) - 1):
            self.update_asteroid()
            self.update_planet()

    def visualize(self):
        """
        Creates a histogram of the period distribution for the asteroids.
        """
        # Orbital period histogram
        dis_sun = np.sqrt(np.sum((self.asteroids_pos - self.sun_pos)**2,  axis=1))[:,None]
        plotlist = np.sqrt(dis_sun**3)/self.const.orbital_period_jup
        plt.figure(1)
        plt.hist(plotlist, edgecolor="black", bins=25)

        # Distance to sun histogram, jupiter's semi major axis included
        plt.figure(2)
        plt.axvline(self.const.smaxis_jup, label="Jupiter SMA", linewidth=2, color='#CC3030')
        plt.hist(dis_sun, edgecolor="black", bins=25)
        plt.legend(fontsize=14, frameon=True, fancybox=True, edgecolor="#000066")
        plt.show()



if __name__ == "__main__":
    c = Constants()
    #sun = Astro_body(c.initial_pos_sun, c.initial_vel_sun, c.MSOL, 0)
    #jupiter = Astro_body(c.initial_pos_jup, c.initial_vel_jup, c.mass_jup, c.ecc_jup)

    #total_time, time_step, amount_of_asteroids)
    test = Kirkwood_solver(100, 0.002, 500, c)
    test.run_N_body_sim()
    print "sun",test.sun_pos
    print "jup",test.jup_pos
    #print "ast",test.asteroids_pos
    print "len", len(test.asteroids_pos)
    test.visualize()

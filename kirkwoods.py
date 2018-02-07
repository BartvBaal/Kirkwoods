"""
    Computational Astrophysics: Numerical Integration final assignment
    Contributors:               Bart van Baal, Boris Wolvers
    Filename:                   kirkwoods.py
    Description:
"""

import numpy as np
import os
import math
import matplotlib.pyplot as plt
import random
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
        self.GM = 4*math.pi*math.pi
        self.MSOL = 1

        # Jupiter info
        self.smaxis_jup = 5.2044
        self.orbital_period_jup = self.smaxis_jup**(3/2.)
        self.mass_jup = 1/1047.
        self.ecc_jup = 0.0489

        # Point 0,0 to be the center of mass of sun/jupiter system
        # Not-too-small timesteps will create spiral-in effect on long runs
        self.offset_cm = self.smaxis_jup / (1+(1/self.mass_jup))

        # Initial velocity of jupiter
        self.start_vel_jup = math.sqrt(
                ((self.GM)/self.smaxis_jup)*((1+self.ecc_jup)/(1-self.ecc_jup)))

        # Initial conditions of jupiter and sun (x0, y0, z0) and (vx0, vy0, vz0)
        self.initial_pos_jup = np.array([0,
                            self.smaxis_jup*(1-self.ecc_jup)-offset_cm, 0])
        self.initial_vel_jup = np.array([self.start_vel_jup, 0, 0])

        self.initial_pos_sun = np.array([0, -self.offset_cm, 0])
        self.initial_vel_sun = np.array(
                    [-self.start_vel_jup*self.mass_jup,0, 0])


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

        # Combining constants to slightly improve numerical speed
        self.gm_time_step = self.const.GM*self.time_step
        self.gm_mass_sol_time_step = self.const.GM*self.const.MSOL*self.time_step
        self.gm_massjup_time_step = self.const.GM*self.const.mass_jup*self.time_step

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

        return np.array(asteroids_pos), np.array(asteroids_vel)


    def update_planet(self):
        """
        Updates the location and velocity for a planet around the sun for a
        single timestep.
        """
        # Determine distance between sun and jupiter
        distance = math.sqrt(np.sum((self.sun_pos - self.jup_pos)**2))

        # Update solar velocity and location; vdot = GM*vector/(distance**3)
        self.sun_vel = (self.sun_vel -
            self.gm_massjup_time_step*(self.sun_pos - self.jup_pos) / distance**3)
        self.sun_pos = self.sun_pos + self.sun_vel*self.time_step

        # Update jupiter velocity and location; vdot = GM*vector/(distance**3)
        self.jup_vel = (self.jup_vel -
            self.gm_mass_sol_time_step*(self.jup_pos - self.sun_pos) / distance**3)
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

        # Check and remove for runaway asteroids
        index_far_asteroids = np.where(dis_sun > 7)[0]
        if len(index_far_asteroids) >= 1:
            dis_sun = np.delete(dis_sun, index_far_asteroids, axis=0)
            dis_jup = np.delete(dis_jup, index_far_asteroids, axis=0)

            self.asteroids_pos = np.delete(self.asteroids_pos, index_far_asteroids, axis=0)
            self.asteroids_vel = np.delete(self.asteroids_vel, index_far_asteroids, axis=0)

        # Update asteroids velocities and locations; vdot = GM*vector/(distance**3)
        self.asteroids_vel = (self.asteroids_vel -
           self.gm_time_step*((self.const.MSOL*(self.asteroids_pos - self.sun_pos) / (dis_sun**3)) +
           (self.const.mass_jup*(self.asteroids_pos - self.jup_pos) / (dis_jup**3))))
        self.asteroids_pos = self.asteroids_pos + self.asteroids_vel*self.time_step


    def run_N_body_sim(self, display=False, start_frac=0.0, update_point=100):
        """
        Calls the update_asteroid and update_planet functions in order to run
        the simulation for the amount of asteroids specified in the init.
        See those specific functions for more comments on their workings.
        If display is set to True, a live feed will be given for the first
        thirty asteroids. start_frac is the point at which the visualization
        starts (so for 0.3 it will start after 30% of the run has been done),
        while update_point is how many steps are in between each update.
        """
        # Save the initial semimajor axis distribution
        self.initial_smaxis_asteroids = (self.find_smaxis_asteroids() /
                                         self.const.smaxis_jup)

        # For live visualization
        if display:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            plt.pause(4)  # Delay for a few seconds

        # Perform n_iterations-1 steps (intialization counts for the first step)
        for i in range(int(self.n_iterations) - 1):
            if i%(500/self.time_step) == 0:
                print i*self.time_step
            self.update_asteroid()
            self.update_planet()

            # Only do a live feed at certain conditions
            if display:
                if i % update_point == 0 and i > start_frac*(int(self.n_iterations)):
                    # Remove the last view so we only see the current asteroid positions
                    plt.cla()
                    ax.scatter(*self.asteroids_pos.T, c="#3399FF", s=1)  # Transpose them
                    ax.scatter(*self.jup_pos.T, c="#660000", s=115)
                    ax.scatter(*self.sun_pos.T, c="#FFA31A", s=250)
                    ax.set_xlim3d(-6, 6)
                    ax.set_xlabel("X (AU)")
                    ax.set_ylim3d(-6, 6)
                    ax.set_ylabel("Y (AU)")
                    ax.set_zlim3d(-.1, .1)
                    ax.set_zlabel("Z (AU)")
                    ax.set_title("Iteration: {}, Asteroids: {}".format(i, len(self.asteroids_pos)))
                    fig.canvas.draw()
                    plt.pause(0.05)
        if display:
            plt.show()


    def find_smaxis_asteroids(self):
        """
        a = -GM/(2E) where E is total energy: E = v**2/2 - GM/r
        Used for the visualization.
        """
        dis_sun = np.sqrt(np.sum((self.asteroids_pos - self.sun_pos)**2,  axis=1))[:,None]
        energy = np.sum(self.asteroids_vel**2, axis=1)[:,None]/2 - self.const.GM/dis_sun

        return -self.const.GM/(2*energy)


    def save_data(self, inlist, namelist):
        """
        Function to save the generated data. Will save some default lists and
        the items inside inlist as separate .npy files. Will use the names in
        namelist for these additional files, so make sure the two lists are the
        same lenght. Recover the data with the np.load() function.
        """
        total_time = int(self.time_step * self.n_iterations)
        subdir = str(self.number_of_asteroids) + "_Asteroids_" + str(total_time) + "_Time_" + str(self.time_step) + "_Stepsize"
        if not os.path.exists("%s/%s" % ("Results", subdir)):
            os.makedirs("%s/%s" % ("Results", subdir))

        # Save all the object values, and then any elements in the inlist
        np.save("%s/%s/initial.npy" % ("Results", subdir), self.initial_smaxis_asteroids)
        np.save("%s/%s/sunpos.npy" % ("Results", subdir), self.sun_pos)
        np.save("%s/%s/sunvel.npy" % ("Results", subdir), self.sun_vel)
        np.save("%s/%s/juppos.npy" % ("Results", subdir), self.jup_pos)
        np.save("%s/%s/jupvel.npy" % ("Results", subdir), self.jup_vel)
        np.save("%s/%s/asteroidspos.npy" % ("Results", subdir), self.asteroids_pos)
        np.save("%s/%s/asteroidsvel.npy" % ("Results", subdir), self.asteroids_vel)
        for i in range(len(inlist)):
            np.save("%s/%s/%s" % ("Results", subdir, str(namelist[i])+".npy"), inlist[i])


    def visualize(self, saving=False):
        """
        Creates a histogram of the period distribution for the asteroids.
        """
        # The initial distribution - clear first to prevent glitches from visualization setups
        plt.figure(1)
        plt.hist(self.initial_smaxis_asteroids, edgecolor="black", bins=np.linspace(0.3, 1., 60))
        plt.title("Initial semi major axis")

        # Orbital period histogram
        smaxis_asteroids = self.find_smaxis_asteroids()
        plotlist = np.sqrt(smaxis_asteroids**3)/self.const.orbital_period_jup

        if saving:
            self.save_data([plotlist], ["plotlist"])  # Do this before plots show/glitch

        plt.figure(2)
        bins = np.linspace(0.3, 1.1, 70)
        plt.hist(plotlist, edgecolor="black", histtype="step", bins=bins, label="Final", lw=2.5, color="#660066", range=(bins.min(),bins.max()))
        plt.hist(self.initial_smaxis_asteroids, histtype="step", bins=np.linspace(0.3, 1.1, 70), label="Initial", lw=2.5, color="#009900")
        plt.legend(fontsize=13, frameon=True, fancybox=True, edgecolor="#000066")

        # Distance to sun histogram, jupiter's semi major axis included
        plt.figure(3)
        plt.axvline(self.const.smaxis_jup, label="Jupiter SMA", linewidth=2, color='#CC3030')
        plt.hist(smaxis_asteroids, edgecolor="black", bins=np.linspace(1.5, 5.5, 150))
        plt.legend(fontsize=14, frameon=True, fancybox=True, edgecolor="#000066")
        plt.show()


if __name__ == "__main__":
    c = Constants()
    
    #total_time, time_step, amount_of_asteroids)
    test = Kirkwood_solver(10000, 1/1024., 25000, c)
    test.run_N_body_sim(display=False)  # Set display to True for live feed
    print "sun",test.sun_pos
    print "jup",test.jup_pos
    #print "ast",test.asteroids_pos
    print "len", len(test.asteroids_pos)
    print c.offset_cm
    test.visualize(saving=True)  # Set saving to True to save the final data

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

# Global parameters - G is PER SOLAR MASS
GM = 4*np.pi**2
MSOL = 3.33e5


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
        
        # Lists of lists of locations & speeds, useful for plotting (for now)
        self.pos = [[self.x0], [self.y0]]
        self.vel = [[self.vx0], [self.vy0]]

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

class Simulation(object):
    """
    Class which executes the simulation, from initialization to visualization
    """
    def __init__(self, amount_of_asteroids, total_time, time_step):
        """
        NOTE: CURRENTLY 2D!
        Always initializes the Sun and Jupiter
        TODO: figure out what decent masses/speeds/locations for asteroids are
              does mass even matter?
        """
        # Jupiter info:
        # semi-major axis: 5.2044 AU
        # eccentricity:    0.0489
        # orbital period:  11.862 year
        # orbital speed:   13.07 km/s => 2.758 AU/year (km/s : AU/y = 1:0.210945021)
        # mass:            1/1047 SolarMass
        
        # Sun at the center, no initial speeds
        self.Sun = Kirkwoods([0, 0], [0, 0], MSOL, total_time, time_step)
        
        jup_a = 5.2044
        jup_m = 1/1047.

        self.Jupiter = Kirkwoods([0, jup_a],
                                 [np.sqrt((GM)/jup_a), 0],
                                 jup_m, total_time, time_step)
        
        ast_sema = (2., 5.)
        ast_mass = 0  # Temp value but it's irrelevant? even for CoM it doesn't really impact stuff

        self.asteroids = []  # List off asteroids, (for now)
        
        # Currently always starts with x=0, vy=0 - thinking of mixing this up (is that needed?)
        for amount in range(amount_of_asteroids):
            startloc = random.uniform(*ast_sema)
            startvel = np.sqrt(GM/startloc)
            asteroid = Kirkwoods([0, startloc],
                                 [startvel, 0],
                                 ast_mass, total_time, time_step)
            self.asteroids.append(asteroid)

    def run_two_body_sim(self, body1, body2):
        """
        Simulates just Sun&Jupiter -> start for 3bodys
        TODO: this should become the first step in the multi_body sims, where
              we first update Jupiter to the next step and then all asteroids
        """
        Mtot = body1.mass + body2.mass
        dimensions = len(body1.pos)

        # Update the location and speeds for these two bodies
        # In case distance doesn't stay the same recalculate inside loop
        # Note, can use this distance to throw away particles i suppose
        
        ## Semi working - doesn't use center of mass & sun doesnt move 
        ## (could simply pretend the sun doesnt actually move, saves calculations...)
        # Distance only checks x&y as Jupiter moves along the z=0 plane
        for i in range(len(body1.time_array)):
            distance = ((body1.pos[0][-1]-body2.pos[0][-1])**2 + 
                        (body1.pos[1][-1]-body2.pos[1][-1])**2)**.5
            for dim in range(dimensions):
                spe1 = body1.vel[dim][-1] - (GM*body1.pos[dim][-1] / 
                                            (distance**3))*body1.time_step
                loc1 = body1.pos[dim][-1] + spe1*body1.time_step
                
                body1.pos[dim].append(loc1)
                body1.vel[dim].append(spe1)

                spe2 = body2.vel[dim][-1] - (GM*body2.pos[dim][-1] /
                                            (distance**3))*body2.time_step
                loc2 = body2.pos[dim][-1] + spe2*body2.time_step

                body2.pos[dim].append(loc2)
                body2.vel[dim].append(spe2)
#            print distance

            # Calculates new pos and velo using Euler-Cromer method
#            vi_plus_1 = velocity[-1] - self.time_step*position[-1]
#            xi_plus_1 = position[-1] + self.time_step*vi_plus_1

#            position.append(xi_plus_1)
#            velocity.append(vi_plus_1)

#        return np.asarray(position), np.asarray(velocity)
    
    
    
    

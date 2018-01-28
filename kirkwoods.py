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
MSOL = 1


class Kirkwoods(object):
    """
    This class contains default initialisations bla bla
    """

    def __init__(self, position, velocity, mass, eccentricity,
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

        # Initial mass & eccentricity
        self.mass = mass
        self.ecc = eccentricity

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
        
        # Point 0,0 to be the center of mass
        # Calculate offset from CoM
        jup_a = 5.2044
        jup_m = 1/1047.
        jup_e = 0.0489
        com_offset = jup_a / (1+(1/jup_m))

        self.Sun = Kirkwoods([0, -com_offset],
                             [0, 0],  # Unsure how to set initial sun speed
                             MSOL, 0,
                             total_time, time_step)

        self.Jupiter = Kirkwoods([0, jup_a*(1-jup_e)-com_offset],
                                 [np.sqrt(((GM)/jup_a)*((1+jup_e)/1-jup_e)), 0],
                                 jup_m, jup_e,
                                 total_time, time_step)
        
        ast_sema = (2., 5.)
        ast_mass = 0  # Temp value but it's irrelevant? even for CoM it doesn't really impact stuff

        self.asteroids = []  # List off asteroids (for now)
        
        # Currently always starts with x=0, vy=0 - thinking of mixing this up (is that needed?)
        for amount in range(amount_of_asteroids):
            startecc = 0  # No eccentricity for now
            startloc = random.uniform(*ast_sema)*(1-startecc)
            startvel = np.sqrt(((1+startecc)/(1-startecc))*GM/startloc)
            asteroid = Kirkwoods([0, startloc],
                                 [startvel, 0],
                                 ast_mass, 0,
                                 total_time, time_step)
            self.asteroids.append(asteroid)

    def update_planet(self, sun, planet, dimensions):
        """
        Works but has some very naive assumptions...
        Updates the location and velocity for a planet around the sun for a
        single timestep.
        Only has a 2D distance as they move in the z=0 plane
        """
        distance = ((sun.pos[0][-1]-planet.pos[0][-1])**2 + 
                    (sun.pos[1][-1]-planet.pos[1][-1])**2)**.5
        
        for dim in range(dimensions):
            sunv = sun.vel[dim][-1] - (GM*sun.pos[dim][-1] / 
                                      (distance**3))*sun.time_step
            sunl = sun.pos[dim][-1] + sunv*sun.time_step
            
            sun.pos[dim].append(sunl)
            sun.vel[dim].append(sunv)

            planetv = planet.vel[dim][-1] - (GM*planet.pos[dim][-1] /
                                            (distance**3))*planet.time_step
            planetl = planet.pos[dim][-1] + planetv*planet.time_step

            planet.pos[dim].append(planetl)
            planet.vel[dim].append(planetv)

    def run_two_body_sim(self, body1, body2):
        """
        Testfunction, currently working on replacement
        Simulates just Sun&Jupiter -> start for 3bodys
        TODO: this should become the first step in the multi_body sims, where
              we first update Jupiter to the next step and then all asteroids
        """
        Mtot = body1.mass + body2.mass
        dimensions = len(body1.pos)

        for i in range(len(body1.time_array)):
            self.update_planet(body1, body2, dimensions)

    def run_three_body_sim(self, body1, body2, body3):
        """
        Runs a three body simulation. Plan to expand this to an N body system
        TODO: resolve the issue with the different masses and what to input...
        """
        Mtot = body1.mass + body2.mass + body3.mass
        dimensions = len(body1.pos)  # Currently still 2D
        time_step = body1.time_step  # Currently all have the same timestep

        # Only 2D for now!!
        for i in range(len(body1.time_array)):
            # Get the distances between the three objects
            sunplanet = ((body1.pos[0][-1]-body2.pos[0][-1])**2 + 
                         (body1.pos[1][-1]-body2.pos[1][-1])**2)**.5
            sunroid = ((body1.pos[0][-1]-body3.pos[0][-1])**2 + 
                       (body1.pos[1][-1]-body3.pos[1][-1])**2)**.5
            planetroid = ((body2.pos[0][-1]-body3.pos[0][-1])**2 + 
                          (body2.pos[1][-1]-body3.pos[1][-1])**2)**.5

            for dim in range(dimensions):
                # Sun only impacted by planet and vice-versa
                # Set current locations to prevent off-by-one actions
                sun_loc = body1.pos[dim][-1]
                pln_loc = body2.pos[dim][-1]
                ast_loc = body3.pos[dim][-1]
                
                sunv = body1.vel[dim][-1] - (GM*body2.mass*(sun_loc - pln_loc) / 
                                            (sunplanet**3))*time_step
                sunl = body1.pos[dim][-1] + sunv*time_step
                
                body1.pos[dim].append(sunl)
                body1.vel[dim].append(sunv)

                plnv = body2.vel[dim][-1] - (GM*body1.mass*(pln_loc - sun_loc) /
                                            (sunplanet**3))*time_step
                plnl = body2.pos[dim][-1] + plnv*time_step

                body2.pos[dim].append(plnl)
                body2.vel[dim].append(plnv)

                roidv = body3.vel[dim][-1] - ((GM*body1.mass*(ast_loc - sun_loc) /
                                             (sunroid**3))*time_step + 
                                             (GM*body2.mass*(ast_loc - pln_loc) /
                                             (planetroid**3))*time_step)
                roidl = body3.pos[dim][-1] + roidv*time_step

                body3.pos[dim].append(roidl)
                body3.vel[dim].append(roidv)




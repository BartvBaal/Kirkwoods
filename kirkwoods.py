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
        self.z0 = position[2]

        # Initial velocity
        self.vx0 = velocity[0]
        self.vy0 = velocity[1]
        self.vz0 = velocity[2]

        # Lists of lists of locations & speeds, useful for plotting (for now)
        self.pos = [[self.x0], [self.y0], [self.z0]]
        self.vel = [[self.vx0], [self.vy0], [self.vz0]]

        # Initial mass & eccentricity
        self.mass = mass
        self.ecc = eccentricity
        
        # Variable to track the distance to the sun and jupiter
        self.soldi = 0
        self.jupdi = 0

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
        NOTE: switched to 3D - appears to be working correctly
        Always initializes the Sun and Jupiter, gives them starting locations
        and speeds in order to keep the center of mass fixed in the origin 0,0.
        Then it initializes as many asteroids as are specified per
        number_of_asteroids, how long the simulation should last in total_time
        and how big the (initial) time_step of the objects should be. Both the
        time_step and total_time are given in years.
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
        jup_startvel = np.sqrt(((GM)/jup_a)*((1+jup_e)/1-jup_e))

        self.Sun = Kirkwoods([0, -com_offset, 0],
                             [-jup_startvel*com_offset/jup_a, 0, 0],
                             MSOL, 0,
                             total_time, time_step)

        self.Jupiter = Kirkwoods([0, jup_a*(1-jup_e)-com_offset, 0],
                                 [jup_startvel, 0, 0],
                                 jup_m, jup_e,
                                 total_time, time_step)
        
        ast_sema = (2, 5.)
        ast_mass = 0  # Temp value but it's irrelevant? even for CoM it doesn't really impact stuff

        self.asteroids = []  # List off asteroids (for now)
        
        # Create all the asteroid Kirkwoods objects
        for amount in range(amount_of_asteroids):
            startecc = 0  # No eccentricity for now
            startloc = random.uniform(*ast_sema)*(1-startecc)
            startvel = np.sqrt(((1+startecc)/(1-startecc))*GM/startloc)

            # Randomize the starting point of the asteroid's orbit
            orb_loc = random.uniform(0, 2*np.pi)
            startx = np.cos(orb_loc)*startloc
            starty = np.sin(orb_loc)*startloc
            startz = 0  # Start in the jupiter-sun plane
            startvelx = np.sin(orb_loc)*startvel
            startvely = -np.cos(orb_loc)*startvel
            startvelz = random.uniform(-startvel, startvel)/100  # Temp range

            # Initialize the asteroid & have set distances to Jupiter and the Sun
            asteroid = Kirkwoods([startx, starty, startz],
                                 [startvelx, startvely, startvelz],
                                 ast_mass, startecc,
                                 total_time, time_step)
            asteroid.sundi = self.get_distance(self.Sun, asteroid)
            asteroid.jupdi = self.get_distance(self.Jupiter, asteroid)
            self.asteroids.append(asteroid)

    def get_distance(self, body1, body2):
        """
        Function to get the distance between two bodies
        """
        dimensions = len(body1.pos)
        distance = 0
        for dim in range(dimensions):
            distance += (body1.pos[dim][-1] - body2.pos[dim][-1])**2
        return np.sqrt(distance)

    def update_planet(self, sun, planet, dimensions):
        """
        Updates the location and velocity for a planet around the sun for a
        single timestep.
        Only has a 2D distance as they move in the z=0 plane
        """
        # Sun and Jupiter should never get a z-location != 0 so only 2D
        distance = ((sun.pos[0][-1]-planet.pos[0][-1])**2 + 
                    (sun.pos[1][-1]-planet.pos[1][-1])**2)**.5

        # Currently does 3D although it only needs to do 2D - might save time to hardcode range(2)
        for dim in range(dimensions):
            sun_loc = sun.pos[dim][-1]
            pln_loc = planet.pos[dim][-1]
            
            sunv = sun.vel[dim][-1] - (GM*planet.mass*(sun_loc - pln_loc) / 
                                      (distance**3))*sun.time_step
            sunl = sun.pos[dim][-1] + sunv*sun.time_step
            
            sun.pos[dim].append(sunl)
            sun.vel[dim].append(sunv)

            planetv = planet.vel[dim][-1] - (GM*sun.mass*(pln_loc - sun_loc) /
                                            (distance**3))*planet.time_step
            planetl = planet.pos[dim][-1] + planetv*planet.time_step

            planet.pos[dim].append(planetl)
            planet.vel[dim].append(planetv)

    def update_asteroid(self, body1, body2, body3, dimensions):
        """
        Updates the position for the asteroid *before* the sun and planet have
        been updated. Body1 as star, body2 as planet, body3 as asteroid.
        dimensions should be the same for all objects. Currently works for 3D.
        """
        for dim in range(dimensions):
                # Set locations, should update asteroids before updating
                # sun/planet to prevent off-by-one errors
                sun_loc = body1.pos[dim][-1]
                pln_loc = body2.pos[dim][-1]
                ast_loc = body3.pos[dim][-1]
                
                roidv = body3.vel[dim][-1] - ((GM*body1.mass*(ast_loc - sun_loc) /
                                              (body3.sundi**3))*body3.time_step + 
                                              (GM*body2.mass*(ast_loc - pln_loc) /
                                               (body3.jupdi**3))*body3.time_step)
                roidl = body3.pos[dim][-1] + roidv*body3.time_step

                body3.pos[dim].append(roidl)
                body3.vel[dim].append(roidv)
        body3.sundi = self.get_distance(body1, body3)
        body3.jupdi = self.get_distance(body2, body3)

    def run_two_body_sim(self, body1, body2):
        """
        DEPRICATED SHOULD NO LONGER BE USED
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
        SHOULD BE REPLACED BY run_N_body_sim FUNCTION COMPLETELY!
        Runs a three body simulation. Plan to expand this to an N body system
        So for an N body simulation it should update the planet once then do
        all the asteroids. Might be easier if asteroids are natively aware how
        far from the sun&planet they are (so the Kirkwoods object should be expanded)
        Once that is added it's no longer needed to calculate sunroid/planetroid
        for each of the asteroids as they already know this & that solves an issue
        of resetting values constantly. So then it becomes update_planet -> update_all_asteroids
        """
        dimensions = len(body1.pos)  # Currently still 2D

        for i in range(len(body1.time_array)):
            self.update_planet(body1, body2, dimensions)
            self.update_asteroid(body1, body2, body3, dimensions)

    def run_N_body_sim(self, body1, body2, body3list):
        """
        Does the three-body simulation for all objects in the body3list, but
        will only update the positions of body1 and body2 once per timestep.
        Note that time_step is build into the objects themselves.
        body1 and body2 are the main gravitational powers in the field, while
        body3list is a list containing all lesser gravitational objects which
        should be updated as a result of the other two objects.
        """
        dimensions = len(body1.pos)  # Currently still 2D
        
        for i in range(len(body1.time_array)):
            for body in body3list:
                self.update_asteroid(body1, body2, body, dimensions)
            self.update_planet(body1, body2, dimensions)

### OLD CODE BELOW ITS PUT INTO FUNCTIONS, ONLY HERE AS TEMPORARY BACKLOG ###
#            for dim in range(dimensions):
#                # Set locations, as planet & sun already got updated it's -2
#                sun_loc = body1.pos[dim][-2]
#                pln_loc = body2.pos[dim][-2]
#                ast_loc = body3.pos[dim][-1]
#                
#                roidv = body3.vel[dim][-1] - ((GM*body1.mass*(ast_loc - sun_loc) /
#                                              (sunroid**3))*time_step + 
#                                              (GM*body2.mass*(ast_loc - pln_loc) /
#                                               (planetroid**3))*time_step)
#                roidl = body3.pos[dim][-1] + roidv*time_step

#                body3.pos[dim].append(roidl)
#                body3.vel[dim].append(roidv)

#        # Only 2D for now!!
#        for i in range(len(body1.time_array)):
#            # Get the distances between the three objects
#            sunplanet = ((body1.pos[0][-1]-body2.pos[0][-1])**2 + 
#                         (body1.pos[1][-1]-body2.pos[1][-1])**2)**.5
#            sunroid = ((body1.pos[0][-1]-body3.pos[0][-1])**2 + 
#                       (body1.pos[1][-1]-body3.pos[1][-1])**2)**.5
#            planetroid = ((body2.pos[0][-1]-body3.pos[0][-1])**2 + 
#                          (body2.pos[1][-1]-body3.pos[1][-1])**2)**.5

#            for dim in range(dimensions):
#                # Sun only impacted by planet and vice-versa
#                # Set current locations to prevent off-by-one actions
#                sun_loc = body1.pos[dim][-1]
#                pln_loc = body2.pos[dim][-1]
#                ast_loc = body3.pos[dim][-1]
#                
#                sunv = body1.vel[dim][-1] - (GM*body2.mass*(sun_loc - pln_loc) / 
#                                            (sunplanet**3))*time_step
#                sunl = body1.pos[dim][-1] + sunv*time_step
#                
#                body1.pos[dim].append(sunl)
#                body1.vel[dim].append(sunv)

#                plnv = body2.vel[dim][-1] - (GM*body1.mass*(pln_loc - sun_loc) /
#                                            (sunplanet**3))*time_step
#                plnl = body2.pos[dim][-1] + plnv*time_step

#                body2.pos[dim].append(plnl)
#                body2.vel[dim].append(plnv)

#                roidv = body3.vel[dim][-1] - ((GM*body1.mass*(ast_loc - sun_loc) /
#                                             (sunroid**3))*time_step + 
#                                             (GM*body2.mass*(ast_loc - pln_loc) /
#                                             (planetroid**3))*time_step)
#                roidl = body3.pos[dim][-1] + roidv*time_step

#                body3.pos[dim].append(roidl)
#                body3.vel[dim].append(roidv)




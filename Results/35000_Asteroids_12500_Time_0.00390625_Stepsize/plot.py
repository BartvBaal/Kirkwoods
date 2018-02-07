import numpy as np
import os
import math
import matplotlib.pyplot as plt
from matplotlib import animation
import random
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# Change border size of plots
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5

initial_smaxis_asteroids = np.load("initial.npy")
asteroids_pos = np.load("asteroidspos.npy")
asteroids_vel = np.load("asteroidsvel.npy")
sun_pos = np.load("sunpos.npy")
sun_vel = np.load("sunvel.npy")
jup_pos = np.load("juppos.npy")
jup_vel = np.load("jupvel.npy")
plotlist = np.load("plotlist.npy")
smaxis_jup = 5.2044
jup_period = 11.862
GM = 4*np.pi*np.pi

dis_sun = np.sqrt(np.sum((asteroids_pos - sun_pos)**2,  axis=1))[:,None]
energy = np.sum(asteroids_vel**2, axis=1)[:,None]/2 - GM/dis_sun

smaxis_asteroids = -GM/(2*energy)
newjupaxis = -GM/(2*(np.sum(jup_vel**2)/2 - GM/np.sum((jup_pos - sun_pos)**2)**.5))
print newjupaxis

index_weird_asteroids = np.where(smaxis_asteroids < 0)[0]
print index_weird_asteroids
print dis_sun[555],  energy[555]
index_close_asteroids = np.where(dis_sun < .3)[0]
print index_close_asteroids

print len(smaxis_asteroids), len(asteroids_pos)
twotoone_res_R = (5.2044**(3/2.) * .5)**(2/3.)
threetoone_res_R = (5.2044**(3/2.) / 3)**(2/3.)
circlebins = np.linspace(0, 2*np.pi, 1000)
print (smaxis_jup**(3/2.)/2)**(2/3.)
first_orbit = [np.cos(circlebins)*twotoone_res_R, np.sin(circlebins)*twotoone_res_R, circlebins*0]
second_orbit = [np.cos(circlebins)*threetoone_res_R, np.sin(circlebins)*threetoone_res_R, circlebins*0]

# The initial distribution - clear first to prevent glitches from visualization setups
plt.figure(1)
plt.hist(initial_smaxis_asteroids*smaxis_jup, edgecolor="black", bins=np.linspace(2., 5., 61))
plt.xlabel("Semi-major Axis (AU)", fontsize=17)
plt.ylabel("Number of Asteroids (per .05 AU bin)", fontsize=17)
plt.xticks(np.linspace(2, 5, 13))
plt.tick_params(labelsize=14, width=1.5)
plt.title("Initial semi major axis", fontsize=20)

plt.figure(2)
bins = np.linspace(0.3, 1.1, 70)
plt.hist(plotlist, edgecolor="black", histtype="step", bins=bins, label="Final", lw=2.5, color="#660066", range=(bins.min(),bins.max()))
plt.hist(initial_smaxis_asteroids, histtype="step", bins=np.linspace(0.3, 1.1, 70), label="Initial", lw=2.5, color="#009900")
plt.legend(fontsize=13, frameon=True, fancybox=True, edgecolor="#000066")

# Distance to sun histogram, jupiter's semi major axis included
plt.figure(3)
plt.hist(smaxis_asteroids, edgecolor="black", bins=np.linspace(1.9, 4.3, 481))
plt.axvline((jup_period/2.)**(2/3.), label="2:1 Resonance", linewidth=2, color='#CC3030', ls="dashed")
plt.axvline((jup_period/3.)**(2/3.), label="3:1 Resonance", linewidth=2, color='#660066', ls="dotted")
plt.axvline((jup_period/2.5)**(2/3.), label="5:2 Resonance", linewidth=2, color='#CC7A00', ls="-.")
plt.legend(fontsize=24, frameon=True, fancybox=True, edgecolor="#000066")
plt.xlim([1.9, 4.3])
plt.text(2.52, 405, "3:1", fontsize=15, color="#660066")
plt.text(3.29, 405, "2:1", fontsize=15, color="#CC3030")
plt.text(2.84, 405, "5:2", fontsize=15, color="#CC7A00")
plt.tick_params(labelsize=17, width=1.5)
plt.xlabel("Semi-major Axis (AU)", fontsize=24)
plt.ylabel("Number of Asteroids (per .005 AU bin)", fontsize=24)
plt.xticks(np.linspace(2, 4.25, 10))
plt.title("Asteroid belt simulation and the Kirkwood gaps, 24977/35000 remaining, 12500 years", fontsize=28)
plt.show()

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#plt.pause(5)
#plt.cla()
#ax.scatter(*asteroids_pos.T, c="#3399FF", s=1)  # Transpose them
#ax.scatter(*jup_pos.T, c="#660000", s=115)
#ax.scatter(*sun_pos.T, c="#FFA31A", s=250)
#ax.plot(*first_orbit, c="Black", lw=2, ls="dashed")
#ax.plot(*second_orbit, c="Black", lw=2, ls="dotted")
#ax.set_xlim3d(-6, 6)
#ax.set_xlabel("X (AU)")
#ax.set_ylim3d(-6, 6)
#ax.set_ylabel("Y (AU)")
#ax.set_zlim3d(-.1, .1)
#ax.set_zlabel("Z (AU)")
#ax.set_title("Iteration: {}, Asteroids: {}".format("End", len(asteroids_pos)))
#fig.canvas.draw()
#plt.show()



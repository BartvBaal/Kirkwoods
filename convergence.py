import numpy as np
import matplotlib.pyplot as plt

convergence = np.load("convergence.npy")
steps = [1/64., 1/128., 1/256., 1/512., 1/1024., 1/2048.]

print np.abs(convergence-5.2044)

plt.loglog(steps, np.abs(convergence-5.2044))
plt.xlabel("Stepsize")
plt.ylabel("Difference from initial Semi-major axis (AU)")
plt.title("Jupiter's decaying orbit for 25000 years")
plt.show()

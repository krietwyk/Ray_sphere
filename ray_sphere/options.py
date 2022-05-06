from ray_sphere.plot_chamber import plot_sphere
import numpy as np
import matplotlib.pyplot as plt

plt.interactive(True)

flag_plot = 1

r0 = 1
P0 = np.array([0, 0, 0])
if flag_plot == 1:
    ax = plot_sphere(r0, P0, 50, 0)
from ray_sphere.plot_chamber import plot_sphere
import numpy as np
import matplotlib.pyplot as plt

from ray_sphere.distances import norm_e
from ray_sphere.photon import * 
# from ray_sphere.enclosure import chamber, port_create, make_port_onaxis, \
#     make_port_offaxis

plt.interactive(True)
flag_plot = 1

# Define the chamber parameters here and have this loaded into each of the 
# relevant modules
# Top modules are enclosure and plot chamber, these require values defined r0,
# P0 and ax, so these must be specifically defined but can be left as global
# variables elsewhere.

######################## Define sphere
r0 = 25
rhos = 0.9
P0 = np.array([0, 0, 0])

if flag_plot == 1:
    fig = plt.figure(0)
    fig.clf()
    ax = plt.axes(projection='3d')
    plot_sphere(r0, P0, 50, ax, 0)

###################### Define ports    
### Port 1
n1 = np.array([0, 0, -1])
p1 = np.array([0, 0, -1])
r1 = 10
h1 = 5

# Port 1 end
rd1 = -1
rho1 = -1
NA1 = 1
ni1 = 1

### Port 2
n2 = np.array([0, 0.2, 0.8])
p2 = np.array([0, 0, 1])
r2 = 10
h2 = 15

# Port 2 end
rd2 = 5
rho2 = 0.9
NA2 = 0.39
ni2 = 1

   
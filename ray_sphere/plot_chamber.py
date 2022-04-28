import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from numpy import pi

def Plot_sphere(r0, P0, npts, figno):
    """Create datapoints for plotting a sphere, t and f are theta and phi.""" 
    f, t = np.mgrid[0:2*np.pi:npts*1j, 0:np.pi:npts*1j]
    x = r0*np.cos(f)*np.sin(t) + P0[0]
    y = r0*np.sin(f)*np.sin(t) + P0[1]
    z = r0*np.cos(t) + P0[2]
    # Create figure for sphere
    fig = plt.figure(figno)
    fig.clf()
    ax = plt.axes(projection='3d')
    # Turn off grid lines
    ax.grid(False)
    # Draw sphere
    ax.plot_surface(x, y, z, color='0.0', alpha=0.3)    
    # Add marker in centre
    ax.scatter(P0[0], P0[1], P0[2]) 
    return ax
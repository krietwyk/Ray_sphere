import matplotlib.pyplot as plt
import numpy as np
from ray_sphere.distances import disp_sph
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

"""
Capitals define points
n and u will always be unit vectors
O and u relate exclusively to photon position and direciton
n is used for the port axis

P0 and r0 are the sphere centre position and radius
p: vector from sphere centre
f and t are used to denote angles phi and theta 
r - cylinder radius
h - cylinder length
P, Q: Points 
v: vector 
"""

def plot_sphere(r0, P0, npts, ax, figno=0):
    """Create datapoints for plotting a sphere, t and f are theta and phi.""" 
    f, t = np.mgrid[0:2*np.pi:npts*1j, 0:np.pi:npts*1j]
    x = r0*np.cos(f)*np.sin(t) + P0[0]
    y = r0*np.sin(f)*np.sin(t) + P0[1]
    z = r0*np.cos(t) + P0[2]
    # Create figure for sphere
    
    # fig = plt.figure(figno)
    # fig.clf()
    # ax = plt.axes(projection='3d')
    # Turn off grid lines
    ax.grid(False)
    # Draw sphere
    ax.plot_surface(x, y, z, color='0.0', alpha=0.3)    
    # Add marker in centre
    ax.scatter(P0[0], P0[1], P0[2]) 
    return #ax

def plot_cylinder(r0, P0, P, n, r, h, npts, figno, ax):
    """Create datapoints for plotting a cylinder intersecting the sphere.
    P - position of port of the cylinder at the sphere
    n - cylinder axis
    r - cylinder radius
    h - cylinder length
    """
    ## replace with create plane?
    # A vector not in the same direction as n
    k = np.array([1, 0, 0])
    if (n == k).all():
        k = np.array([0, 1, 0])
    # A unit vector n1 perpendicular to n
    w1 = np.cross(n, k)
    w1 = w1/np.dot(w1, w1)**0.5
    # A unit vector perpendicular to v and n1
    w2 = np.cross(n, w1)
    w2 = w2/np.dot(w2,w2)
    # Surface ranges over t from 0 to length of axis and f = 0 to 2*pi
    f, d = np.mgrid[0:2*np.pi:npts*1j, 0:h:3j]
    # Generate coordinates for the surface
    X, Y, Z = [P[i] + n[i]*d + r*np.sin(f)*w1[i] + r*np.cos(f)*w2[i] \
               for i in [0, 1, 2]]
    ax.plot_surface(X, Y, Z, color='0', alpha=0.3)
    return

def plot_cylinder_offaxis(r0, P0, p, n, r, h, w1, w2, npts, figno, ax): 
    # Example of intersection between cylinder projection and sphere, Q
    # (v1, disp1, v2, disp2) = disp_sph(w1*r + r0*p-P0, n, P0, r0)
  
    # Plot the intersection between cylinder projection and sphere
    psi_ = np.linspace(0,2*np.pi,npts)
    G = np.zeros([3,npts])
    G2 = np.zeros([3,npts])
    for i in range(len(psi_)):
        # Define a vector for the cylinder end then project towards the sphere, G
        Q = r*(w1*np.cos(psi_[i]) + w2*np.sin(psi_[i])) + r0*p-P0
        (Q1, *_) = disp_sph(Q, n, P0, r0)
        
        # Sphere end of the cylinder
        G[:,i] = Q1
        
        # Far end of the cylinder
        G2[:,i] = Q + n*h
        

        if i > 0:
            x = [G[0, i-1], G2[0, i-1], G2[0,i], G[0,i]]
            y = [G[1, i-1], G2[1, i-1], G2[1,i], G[1,i]]
            z = [G[2, i-1], G2[2, i-1], G2[2,i], G[2,i]]
            verts = [list(zip(x, y, z))]
            ax.add_collection3d(Poly3DCollection(verts))
    ax.plot3D(G[0], G[1], G[2], color = 'black')
    ax.plot3D(G2[0], G2[1], G2[2], color = 'black')
    
    ax.scatter(r0*p[0], r0*p[1], r0*p[2]) # centre of port on the sphere
    

    




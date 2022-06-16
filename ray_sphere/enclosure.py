# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
# import ray_sphere.options

from ray_sphere.distances import norm_e, disp_sph, line_line_int, \
    disp_cylinder, plane_create, disp_line_point
from ray_sphere.plot_chamber import plot_sphere, plot_cylinder, plot_cylinder_offaxis
from ray_sphere.options import *
from ray_sphere.locations import lamb_event
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class chamber:        
    ele = []    
    def add_port(self, port):
        print(f'port number {len(self.ele)} added')
        self.ele.append(port)
 
    def port_empty(self):
        chamber.ele = []

class port_create:
    def __init__(self, n, p, r, h, rd=-1, rho=-1, NA=1, ni=1):
        self.n = n/norm_e(n) # ensure unit vectors
        self.p = p/norm_e(p) # ensure unit vectors
        self.r = r
        self.h = h
        self.P = []
        self.w1 = []
        self.w2 = []              
        self.loc = len(chamber.ele) + 1
        self.rho = rho
        self.rd = rd
        self.NA = NA
        self.ni = ni
        
        # Add the shape that intersects the sphere
        if np.dot(self.n, self.p) == 1: # if on axis
            print('on axis')
            self.port_type = 'on_axis'
            
            (self.P, self.sc) = make_port_onaxis(r0, P0, self.r, self.n)
            self.w1, self.w2 = plane_create(self.n)

            # Plot cirlce on axis plot_cylinder(P, n, r, h, npts, figno, ax):
            plot_cylinder(r0, P0, self.sc*self.n, self.n, self.r, self.h, 50, 1, ax)        
        else: # off axis
            self.port_type = 'off_axis'
            self.w1, self.w2, self.P = make_port_offaxis(r0, self.p, self.n)    
            print('make off axis')
               
            plot_cylinder_offaxis(r0, P0, self.p, self.n, self.r, self.h, 
                                  self.w1, self.w2, 50, 1, ax)        
        if h > 0:
            self.port_type = self.port_type + '_cyl'
            print('Also make cylinder element')
        
        chamber().add_port(port=self)
        
def make_port_onaxis(r0, P0, r, n):
    """"Create ports on axis to the sphere with vector normals
    should point away from sphere centre
    r0 = sphere_radius 
    P0 = sphere_centre 
    r = port radius
    n = normal of port face closest to sphere, pointing outwards   
    """
    # Normalises the vector
    n = n/norm_e(n)
    
    # distance from centre of sphere to centre of circular port
    sc = (r0**2-r**2)**0.5
    # Centre of circle defining port
    Q = P0 + sc*n
    return (Q, sc)

def make_port_offaxis(r0, p0, n):    
    """Make cylinder on sphere that is off-axis."""
    p0 = p0/norm_e(p0)

    # Direction of port normal
    n = n/norm_e(n)
    if np.dot(n, p0) == 1: # check if on axis
        print("not off axis!")
        return
    # unit vector of u+p, points vaguely toward minor axis of ellipse    
    w = (n+p0)/np.dot(n+p0, n+p0)**0.5    
     
    # w1, w2 = Plane_create(w)
    # ##### Assme that w.u ~= 1, or else on axis 
    w2 = np.cross(w, n)
    w2 = w2/np.dot(w2, w2)**0.5
    w1 = np.cross(w2, n)
    
    if np.dot(w1,w) < 0:
        w1 = np.cross(-w2, n)
    w1 = w1/np.dot(w1, w1)**0.5
    Q = r0*p0
    return (w1, w2, Q)

def replot_all(chamber, r0, P0, npts, ax, figno=0):
    fig = plt.figure(figno)
    fig.clf()
    ax = plt.axes(projection='3d')
    plot_sphere(r0, P0, npts, ax, figno)
    for element in chamber.ele:
        if (element.n == element.p).all():
            plot_cylinder(r0, P0, element.P, element.n, element.r, element.h, 
                          npts, figno, ax)
        else:
            plot_cylinder_offaxis(r0, P0, element.p, element.n, element.r, 
                                  element.h, element.w1, element.w2, 
                                  npts, figno, ax)
    return ax

    
    
    

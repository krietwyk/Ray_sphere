# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from ray_sphere.distances import norm_e, dist_sph, line_line_int, dist_line, \
    dist_cylinder, plane_create
from ray_sphere.plot_chamber import plot_sphere, plot_cylinder, plot_cylinder_offaxis

# plt.interactive(True)

class chamber:
    ele = []
    
    def add_port(self, port):
        print(f'port number {len(self.ele)} added')
        self.ele.append(port)

class port_create:
    # if h = 0, ignore second end, only focus 
    def __init__(self, n, p0, r, h, end_type1=1, end_type2=0):
        self.n = n
        self.p0 = p0
        self.r = r
        self.h = h
        self.P = []
        self.w1 = []
        self.w2 = []              
        self.loc = len(chamber.ele) + 1
        self.end_type1 = end_type1
        
        # Add the shape that intersects the sphere
        if np.dot(n, p0) == 1: # if on axis
            print('on axis')
            self.port_type = 'on_axis'
            
            (self.P, self.sc) = make_port_onaxis(r0, P0, self.r, self.n)
            self.w1, self.w2 = plane_create(self.n)
            
            if flag_plot > 0:
                # Plot cirlce on axis plot_cylinder(P, n, r, h, npts, figno, ax):
                plot_cylinder(self.sc*self.n, self.n, self.r, self.h, 50, 1, ax)        
        else: # off axis
            self.port_type = 'off_axis'
            self.w1, self.w2, self.P = make_port_offaxis(r0, self.p0, self.n)    
            print('make off axis')
            
            if flag_plot > 0:   
                plot_cylinder_offaxis(r0, P0, self.p0, self.n, self.r, self.h, 
                                      self.w1, self.w2, 50, 1, ax)        
        if h > 0:
            self.port_type = self.port_type + '_cyl'
            self.end_type2 = end_type2
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
    # Normalises the vector, i.e. P is a unit vector
    if np.dot(n, n.T)**0.5 != 1:    
        n = n/np.dot(n, n)**0.5 # declare P outside of function
    # distance from centre of sphere to centre of circular port
    sc = (r0**2-r**2)**0.5
    # Centre of circle defining port
    Q = P0 + sc*n
    return (Q, sc)

def make_port_offaxis(r0, p0, n):    
    """Make cylinder on sphere that is off-axis."""
    p0 = p0/np.dot(p0, p0)**0.5

    # Direction of port normal
    n = n/np.dot(n, n)**0.5
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


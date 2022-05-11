# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from ray_sphere.distances import norm_e, dist_sph, line_line_int, dist_line, \
    dist_cylinder, plane_create, disp_line_point
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

def is_in_on(photon, chamber, i):
    
    # Closer to the cylinder axis than the radius
    dist_r, *_ = disp_line_point(chamber.ele[i].P, chamber.ele[i].n,  
                                 photon.O)
    # dist_r, _ = dist_line(photon.O, 
    #                       chamber.ele[i].P, 
    #                       chamber.ele[i].n, 
    #                       photon.u)    
    # Make a point that is the centre of the far end cap of the cylinder
    # Q = chamber.ele[i].sc + chamber.ele[i].n*chamber.ele[i].h
    
    # Displacement of project onto cylinder axis from sphere centre
    a = np.dot(photon.O-P0, chamber.ele[i].p)
    
    # Further than the spherical cap but is it not too far away
    # dist_s = np.dot(photon.O-P0, photon.O-P0)**0.5
    # sc = chamber.ele[i].sc
    # sc = (r0**2-chamber.ele[i].r**2)**0.5
    # Closer than the end cap
    dist_c = chamber.ele[i].sc + chamber.ele[i].h
    if dist_r <= chamber.ele[i].r and a >= chamber.ele[i].sc and a <= dist_c:
        print(f'Within on-axis cylinder {i}')
        is_in = 1
    else:
        print(f'Outside of on-axis cylinder {i}')
        is_in = -99
    return is_in


def is_in_off(photon, chamber, i):
    # flag_cont = 1 # if 0, reduce calcs and run to return faster
    w1 = chamber.ele[i].w1
    w2 = chamber.ele[i].w2
    
    # Photon with the cylinder radius from cylinder axis
    dist_r, *_ = disp_line_point(chamber.ele[i].P, chamber.ele[i].n,  
                                 photon.O)
    
    # Make a point that is the centre of the far end cap of the cylinder
    Q = chamber.ele[i].P + chamber.ele[i].n*chamber.ele[i].h
    
    # Distance of project onto cylinder axis, closer than the end cap
    a = np.dot(photon.O-Q, -chamber.ele[i].n)
    
    # Distance from sphere centre onto port axis
    b = np.dot(photon.O-P0, chamber.ele[i].p)
    
    if norm_e(dist_r) <= chamber.ele[i].r and a >= 0 and b >= r0: 
        print(f'Within off-axis cylinder {i}')
        is_in = 1
    else:
        print(f'Outside of off-axis cylinder {i}')
        is_in = -99
    return is_in

def locator(photon, chamber):
    flag_chk_S = 1
    flag_in_in = -99
    flag_in_off = -99
    for i in range(len(chamber.ele)):
        # is in on-axis port
        if chamber.ele[i].port_type[0:2] == 'on':
            flag_in_in = is_in_on(photon, chamber, i)
            if flag_in_in < 0:
                print(f'not in elemment {i}')
        # Is in off
        elif chamber.ele[i].port_type[0:2] == 'of':
            flag_in_off = is_in_off(photon, chamber, i)
            if flag_in_off < 0:
                print(f'not in elemment {i}')
        # If in port don't check if in sphere 
        if flag_in_in > 0 or flag_in_off > 0:
            print(flag_in_in, flag_in_off)
            print(f'In element {i}')
            flag_chk_S = 0
            break        

    # Is in sphere, this is last because it can be a port and sphere sim.
    if flag_chk_S == 1:
        if norm_e(photon.O - P0) <= r0:
            i = -1 # in sphere
            print(f'in sphere {i}')
        else:
            i = -2 # outside of chamber
            print(f'Outside of chamber {i}')
    return i

import numpy as np

from ray_sphere.options import *

from ray_sphere.distances import *

from ray_sphere.enclosure import *

from ray_sphere.plot_chamber import *

from ray_sphere.photon import *

from ray_sphere.locations import *

                      
def event(photon, chamber, i, comp_ind, Q):
    # Interaction with cylinder wall, lambertian reflelction
    if chamber.ele[i].h > 0 and comp_ind == 1:
        a = chamber.ele[photon.loc]
        b = np.dot(Q - (a.P+a.h*a.n), -a.n)
        lamb_event(photon, comp_ind, Q-(a.P+a.h*a.n-b*a.n), Q, rhos, 1, color='black')

    # Interaction with port end (away from sphere)        
    # elif chamber.ele[i].h > 0 and (comp_ind == 2:
    elif comp_ind == 0 or comp_ind == 2:    
        # negative rho means it pass through through, if moving outwards
        dist_c = norm_e(Q - (chamber.ele[i].h*chamber.ele[i].n+chamber.ele[i].P))
        
        if (comp_ind == 2 and chamber.ele[i].rho < 0):  
            if np.dot(photon.u, chamber.ele[i].n) > 0:
                photon.loc, photon.O = -2, Q 
            elif np.dot(photon.u, chamber.ele[i].n) <= 0:
                photon.loc, photon.O = i, Q
        elif comp_ind == 0:
            if np.dot(photon.u, -chamber.ele[i].n) > 0:
                photon.loc, photon.O = -1, Q
        # if it interacts with end port check it is within the detector radius
        elif dist_c <= chamber.ele[i].rd:
            photon = detect(photon, chamber, comp_ind, Q)
        else: 
            photon = lamb_event(photon, comp_ind, chamber.ele[i].n, 
                                              Q, chamber.ele[i].rho, 1, 
                                              color='black')
            
def detect(photon, chamber, element, Q):
    # NA = ni.asin(theta), ni is refractive index
    theta = np.arccos(np.dot(photon.u, chamber.ele[element].n))
    acceptance_angle = np.arcsin(chamber.ele[element].NA/chamber.ele[element].ni)
    if theta <= acceptance_angle:
        print('photon detected')
        ax.plot3D([Q[0], photon.O[0]], 
                  [Q[1], photon.O[1]], 
                  [Q[2], photon.O[2]], color = 'red')
        photon.state = 1
        # photon.O = Q
        
        # photon.fate_rec(photon, element, 1, theta)
    else:
        photon.state = -2
        # photon.O = Q
        # photon.fate_rec(photon, element, -2, theta)
    photon.O = Q
    photon.fate_rec(photon, element, photon.state, theta)
    return photon

def lamb_event(photon, comp, n, Q, rho, flag_plot, color='black'):
    # flag_new_loc = 0
    # if rho < 0: # pass through (open port)
    #     photon.old.append(photon.O)
    #     photon.O = Q
    #     flag_new_loc = 1
    if np.random.rand() <= rho:
        v, k = lamb_ref(photon.u, n)
        # photon.old.append(photon.O)
        if flag_plot >= 1:
            ax.plot3D([Q[0], photon.O[0]], 
                      [Q[1], photon.O[1]], 
                      [Q[2], photon.O[2]], color = 'black')
        photon.O, photon.u  = Q, v
    else:
        photon.state = 0
        photon.fate_rec(photon, comp, photon.state, np.nan)
    return photon #, flag_new_loc

def lamb_ref(u, n): 
    # n must be a unit vector but u doesn't need to be.
    # for u and n pointing in the same directions. 
    # Resulting vector points opposite to surface normal, -n. Returns unit vector.
    
    t = np.pi*np.random.rand() #% acos(1-2*rand);
    f = 2*np.pi*np.random.rand() #% acos(1-2*rand);
    st = np.sin(t)
    
    # k = np.dot(photon.u, n)*n - photon.u
    k = -n
    k = k/norm_e(k)
    v = [k[0]+st*np.cos(f), k[1]+st*np.sin(f), k[2]+np.cos(t)]
    v = v/norm_e(v)
    return v, k

def spec_ref(u, n): 
    return 2*np.dot(n, u)*n - u

        
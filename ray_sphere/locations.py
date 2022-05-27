import numpy as np
from ray_sphere.distances import norm_e, disp_sph, line_line_int, \
    disp_cylinder, plane_create, disp_line_point, line_plane_int

from ray_sphere.plot_chamber import plot_sphere, plot_cylinder, \
    plot_cylinder_offaxis
    
# from ray_sphere.options import *
from ray_sphere.options import P0, r0, rhos, ax

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
    # w1 = chamber.ele[i].w1
    # w2 = chamber.ele[i].w2
    
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

def in_sph(photon, chamber): # Det. intersection
    while  photon.state == 1 and photon.loc == -1:
        if len(chamber.ele) == 0:
            print('hits wall')
            Q, disp_s, _, _ = disp_sph(photon.O, photon.u, P0, r0)
            if np.random.rand() <= rhos:
                v, k = lamb_ref(photon.u, Q-P0)
                photon.old.append(photon.O)
    
                ax.plot3D([Q[0], photon.O[0]], 
                          [Q[1], photon.O[1]], 
                          [Q[2], photon.O[2]], color = 'black')
                photon.O = Q
                photon.u = v
            else:
                photon.state = 0
        
        for i in range(len(chamber.ele)):
            if np.dot(chamber.ele[i].n, chamber.ele[i].p) == 1:
                disp, Q = line_plane_int(photon.O, photon.u, chamber.ele[i].P, 
                                         chamber.ele[i].n)
                if disp > 0 and norm_e(Q-chamber.ele[i].P) < chamber.ele[i].r:
                    print('run event')
                    photon.loc = i
                    photon.old.append(photon.O)
                    ax.plot3D([Q[0], photon.O[0]], 
                              [Q[1], photon.O[1]], 
                              [Q[2], photon.O[2]], color = 'black')
                    photon.O = Q
                    # event(photon, chamber, i, comp_ind, Q)
                    break
                else:
                    print('hits wall')
                    Q, disp_s, _, _ = disp_sph(photon.O, photon.u, P0, r0)
                    if np.random.rand() <= rhos:
                        v, k = lamb_ref(photon.u, Q-P0)
                        photon.old.append(photon.O)
                        ax.plot3D([Q[0], photon.O[0]], 
                                  [Q[1], photon.O[1]], 
                                  [Q[2], photon.O[2]], color = 'black')
                        photon.O = Q
                        photon.u = v
                    else:
                        photon.state = 0
                
            elif np.dot(chamber.ele[i].n, chamber.ele[i].p) != 1:
                # print('off axis')
                (Q, disp1, Q2, disp2) = disp_sph(photon.O, photon.u, P0, r0)
                disp_r, *_ = disp_line_point(chamber.ele[i].P, chamber.ele[i].n, Q)
                if norm_e(disp_r) < chamber.ele[i].r:
                    # event(photon, chamber, i, comp_ind, Q)
                    photon.loc = i
                    break
                else: # no interactions with ports, photon hits wall
                    print('hits wall')
                    Q, disp_s, _, _ = disp_sph(photon.O, photon.u, P0, r0)
                    if np.random.rand() <= rhos:
                        v, k = lamb_ref(photon.u, Q-P0)
                        photon.old.append(photon.O)
    
                        ax.plot3D([Q[0], photon.O[0]], 
                                  [Q[1], photon.O[1]], 
                                  [Q[2], photon.O[2]], color = 'black')
                        photon.O = Q
                        photon.u = v
                    else:
                        photon.state = 0
                        print('photon died')
                        
                        
def event(photon, chamber, i, comp_ind, Q):
    
    # if chamber.ele[i].h == 0: # if there is no length, so this is the end cap
    #     if chamber.ele[i].end_type1 == 0:
    #         ax.plot3D([Q[0], photon.O[0]], 
    #                   [Q[1], photon.O[1]], 
    #                   [Q[2], photon.O[2]], color = 'red')
    #         photon.state = 0 
    #         print('photon detected')

    # elif chamber.ele[i].h > 0 and comp_ind == 0:
    #     photon.old.append(photon.O)
    #     ax.plot3D([Q[0], photon.O[0]], 
    #               [Q[1], photon.O[1]], 
    #               [Q[2], photon.O[2]], color = 'black')
    #     photon.O = Q
    #     photon.loc = i
    
    # Interaction with cylinder wall, lambertian reflelction
    if chamber.ele[i].h > 0 and comp_ind == 1:
        a = chamber.ele[photon.loc]
        b = np.dot(Q - (a.P+a.h*a.n), -a.n)
        lamb_event(photon, Q-(a.P+a.h*a.n-b*a.n), Q, rhos, 1, color='black')
        # if np.random.rand() <= rhos:
        #     a = chamber.ele[photon.loc]
        #     b = np.dot(Q - (a.P+a.h*a.n), -a.n)
        #     v, k = lamb_ref(photon.u, Q - (a.P+a.h*a.n-b*a.n))
        #     photon.old.append(photon.O)

        #     ax.plot3D([Q[0], photon.O[0]], 
        #               [Q[1], photon.O[1]], 
        #               [Q[2], photon.O[2]], color = 'black')
        #     photon.O = Q
        #     photon.u = v
        # else:
        #     photon.state = 0 
    
    # Interaction with port end (away from sphere)        
    # elif chamber.ele[i].h > 0 and (comp_ind == 2:
    elif comp_ind == (0 or 2):    
        # negative rho means it pass through through, if moving outwards
        dist_c = norm_e(Q - chamber.ele[i].h*chamber.ele[i].n+chamber.ele[i].P)
        if chamber.ele[i].rho < 0:  
            if np.dot(photon.u, chamber.ele[i].n) > 0:
                photon.loc = -2 
            elif np.dot(photon.u, chamber.ele[i].n) < 0:
                photon.loc = i
        # if it interacts with end port check it is within the detector radius
        elif dist_c <= chamber.ele[i].rd:
            
            photon = detect(photon, chamber, i, Q)
        else: 
            photon = lamb_event(photon, chamber.ele[i].n, 
                                              Q, chamber.ele[i].rho, 1, 
                                              color='black')
            
def detect(photon, chamber, element, Q):
    # NA = ni.asin(theta), ni is refractive index
    theta = np.acos(np.dot(photon.u, chamber.ele[element].n))
    acceptance_angle = np.asin(chamber.ele[element].NA/chamber.ele[element].ni)
    if theta <= acceptance_angle:
        print('photon detected')
        ax.plot3D([Q[0], photon.O[0]], 
                  [Q[1], photon.O[1]], 
                  [Q[2], photon.O[2]], color = 'red')
        photon.O = Q
        photon.state = -1
        # photon.fate_rec()
    return photon

def lamb_event(photon, n, rho, Q, flag_plot, color='black'):
    # flag_new_loc = 0
    # if rho < 0: # pass through (open port)
    #     photon.old.append(photon.O)
    #     photon.O = Q
    #     flag_new_loc = 1
    if np.random.rand() <= rho:
        v, k = lamb_ref(photon.u, n)
        photon.old.append(photon.O)
        if flag_plot >= 1:
            ax.plot3D([Q[0], photon.O[0]], 
                      [Q[1], photon.O[1]], 
                      [Q[2], photon.O[2]], color = 'black')
        photon.O = Q
        photon.u = v
    else:
        photon.state = 0
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


import numpy as np

from ray_sphere.options import *

from ray_sphere.distances import *

from ray_sphere.enclosure import *

from ray_sphere.plot_chamber import *

from ray_sphere.photon import *

from ray_sphere.locations import *

def in_on_cyl(photon, chamber, i):
    """
    Calculate where adn how the photon interacts by in on-axis cylinder by
    calculating the expected intersection of the ends and walls of the 
    cylinder. The photon interacts with the component with the smallest +'ve 
    displacement.
    """
    while photon.state == 1 and photon.loc == i:
        # 1. Intersecting sphere port?
        if np.dot(photon.u, -chamber.ele[i].n) > 0:
            disp_s, A = line_plane_int(photon.O, photon.u, chamber.ele[i].P, 
                                   -chamber.ele[i].n)
        else:
            disp_s, A = -1, -1                
        # 2. Intersecting and cylinder wall?
        (disp_w , B, d, v2) = disp_cylinder(photon.O, photon.u, chamber.ele[i].P,
                                       chamber.ele[i].n, chamber.ele[i].r)
        comp = [disp_s, disp_w]
        comp_pt = [A, B]
        # 3. Intersecting an on-axis port - End cap?
        if chamber.ele[i].h > 0:
            if np.dot(photon.u, chamber.ele[i].n) > 0:
                disp_c, C = line_plane_int(photon.O, photon.u, 
                                        chamber.ele[i].P+chamber.ele[i].h*chamber.ele[i].n,
                                        chamber.ele[i].n)
            else:
                disp_c, C = -1, -1
            comp.append(disp_c)
            comp_pt.append(C)
            
        comp_ind = comp.index(min(i for i in comp if i > 0))
        Q = comp_pt[comp_ind]
        
        event(photon, chamber, i, comp_ind, Q)

# Not correct at all
def in_off_cyl(photon, chamber, i): # Det. intersection
    """
    Calculate where adn how the photon interacts by in off-axis cylinder by
    calculating the expected intersection of the ends and walls of the 
    cylinder. The photon interacts with the component with the smallest +'ve 
    displacement.
    """
    
    while photon.state == 1 and photon.loc == i:
        (a, A, b, B) = disp_cylinder(photon.O, photon.u, chamber.ele[i].P,
                                     chamber.ele[i].n, chamber.ele[i].r)

        # Is norm_e(A-P0) <= r0 and np.dot(A-P0, chamber.ele[i].n) >= 0
        # if so, photon enters the sphere. Calculate sphere entry point
        G = chamber.ele[i].P+chamber.ele[i].h*chamber.ele[i].n
        if norm_e(A-P0) <= r0 and np.dot(A-P0, chamber.ele[i].n) >= 0:
            # Photon enters the sphere. Calculate sphere entry point:
            (Q1, disp1, Q2, disp2) = disp_sph(photon.O, photon.u, P0, r0)
            comp_pt = Q1
            comp_ind = 0
            
        elif np.dot(A - G, -chamber.ele[i].n) > 0:
            # photon hits wall between sphere and end cap.
            comp_pt = A
            comp_ind = 1
            
        elif chamber.ele[i].h > 0:
            # else photon hits end cap
            disp_c, Q = line_plane_int(photon.O, photon.u, 
                                    chamber.ele[i].P+chamber.ele[i].h*chamber.ele[i].n,
                                    chamber.ele[i].n)
            comp_pt = Q 
            comp_ind = 2

        event(photon, chamber, i, comp_ind, comp_pt)

    
# Working, determines what is hit if anything at all but events not covered
def in_out(photon, chamber):
    
    # Does it intercept with any cylinders, outer port ends or the sphere 
    # sphere
    while photon.state == 1 and photon.loc == -2:
        # Need second set of pos. and disp. value 
        dist_list = []
        ele_list = []
        pt_list = []
        comp_list = []
        
        (_, _, Q, disp_s) = disp_sph(photon.O, photon.u, P0, r0)
        
        dist_list.append(disp_s)
        ele_list.append(-1)
        comp_list.append(0)
        pt_list.append(Q)
        for i in range(len(chamber.ele)):    
            # Consider cylinder wall
            if chamber.ele[i].h > 0: # has some cylinder wall    
                # hit cylinder wall?
                (_, _, disp_w, Q) = disp_cylinder(photon.O, photon.u, chamber.ele[i].P,
                                           chamber.ele[i].n, chamber.ele[i].r)
                # position of centre of outer end cap
                G = chamber.ele[i].P+chamber.ele[i].h*chamber.ele[i].n
                # if np.dot(chamber.ele[i].n,  chamber.ele[i].p) == 1:
                    # on axis
                    # cond1: Vector from cyl wall intercept to outer end cap point
                    # towards the sphere, close to the sphere than end cap
                    # cond2: vector from cyl wall intercept to sphere point away 
                    # from sphere centre, i.e, is Q outside of sphere
                if np.dot(Q-G, -chamber.ele[i].n) > 0 and \
                    np.dot(Q-chamber.ele[i].P, chamber.ele[i].n) > 0:
                        # hits cylinder wall
                        dist_list.append(disp_w)
                        ele_list.append(i)
                        comp_list.append(1)
                        pt_list.append(Q)

            # Consider outer port end
            disp_p, Q = line_plane_int(photon.O, photon.u, 
                                       chamber.ele[i].P+chamber.ele[i].h*chamber.ele[i].n, 
                                     chamber.ele[i].n)
            dist_c = norm_e(Q-(chamber.ele[i].P+chamber.ele[i].h*chamber.ele[i].n))
            if dist_c <= chamber.ele[i].r and np.dot(-chamber.ele[i].n, photon.u) > 0:
                # If the photon hits an out port moving against the normal i.e. towards
                dist_list.append(disp_p)
                ele_list.append(i) 
                pt_list.append(Q)
                comp_list.append(2)
                    
        # Figure out which of the components that are intercepted would be hit 
        # first 
        dist_list_pos = [i for i in dist_list if i > 0]
        if len(dist_list_pos) == 0:
            photon.state = 0 
        else:
            comp_ind = dist_list.index(min(dist_list_pos))
            comp_pt = pt_list[comp_ind]
            element = ele_list[comp_ind]
            # Q = comp_pt
            
            if comp_list[comp_ind] == 2 and chamber.ele[element].rho < 0: 
                photon.loc, photon.O = element, comp_pt
            else:
                photon.state = 0 
        
def in_sph(photon, chamber): # Det. intersection
    while  photon.state == 1 and photon.loc == -1:
        if len(chamber.ele) == 0:
            print('hits wall')
            Q, disp_s, _, _ = disp_sph(photon.O, photon.u, P0, r0)
            lamb_event(photon, 1, Q-P0, Q, rhos, 1, color='black')
            
        for i in range(len(chamber.ele)):
            if np.dot(chamber.ele[i].n, chamber.ele[i].p) == 1: #on axis
                disp, Q = line_plane_int(photon.O, photon.u, chamber.ele[i].P, 
                                          chamber.ele[i].n)
                if disp > 0 and norm_e(Q-chamber.ele[i].P) < chamber.ele[i].r:
                    print('run event')
                    photon.loc = i
                    ax.plot3D([Q[0], photon.O[0]], 
                              [Q[1], photon.O[1]], 
                              [Q[2], photon.O[2]], color = 'black')
                    photon.O = Q
                    break
                else:
                    print('hits wall')
                    Q, disp_s, _, _ = disp_sph(photon.O, photon.u, P0, r0)
                    lamb_event(photon, 1, Q-P0, Q, rhos, 1, color='black')

            elif np.dot(chamber.ele[i].n, chamber.ele[i].p) != 1: # off axis
                # print('off axis')
                (Q, disp1, Q2, disp2) = disp_sph(photon.O, photon.u, P0, r0)
                disp_r, *_ = disp_line_point(chamber.ele[i].P, chamber.ele[i].n, Q)
                if norm_e(disp_r) < chamber.ele[i].r:
                    print('run event')
                    photon.loc = i
                    ax.plot3D([Q[0], photon.O[0]], 
                              [Q[1], photon.O[1]], 
                              [Q[2], photon.O[2]], color = 'black')
                    photon.O = Q
                    break
                else: # Doesn't interact with port openings so must hit a wall
                    print('hits wall')
                    Q, disp_s, _, _ = disp_sph(photon.O, photon.u, P0, r0)
                    lamb_event(photon, 1, Q-P0, Q, rhos, 1, color='black') 

        
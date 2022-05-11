import numpy as np
from ray_sphere.distances import norm_e, disp_sph, line_line_int, \
    disp_cylinder, plane_create, disp_line_point
from ray_sphere.plot_chamber import plot_sphere, plot_cylinder, \
    plot_cylinder_offaxis

def event(photon, chamber, i, comp_ind, Q):
    if chamber.ele[i].h == 0: # if there is no length, so this is the end cap
        if chamber.ele[i].end_type1 == 0:
            ax.plot3D([Q[0], photon.O[0]], 
                      [Q[1], photon.O[1]], 
                      [Q[2], photon.O[2]], color = 'red')
            photon.state = 0 
            print('photon detected')

    elif chamber.ele[i].h > 0 and comp_ind == 0:
        photon.old.append(photon.O)
        ax.plot3D([Q[0], photon.O[0]], 
                  [Q[1], photon.O[1]], 
                  [Q[2], photon.O[2]], color = 'black')
        photon.O = Q
        photon.loc = i
            
    elif chamber.ele[i].h > 0 and comp_ind == 1:
        if np.random.rand() <= rhos:
            a = chamber.ele[photon.loc]
            b = np.dot(Q - (a.P+a.h*a.n), -a.n)
            v, k = lamb_ref(ph.u, Q - (a.P+a.h*a.n-b*a.n))
            photon.old.append(photon.O)

            ax.plot3D([Q[0], photon.O[0]], 
                      [Q[1], photon.O[1]], 
                      [Q[2], photon.O[2]], color = 'black')
            photon.O = Q
            photon.u = v
        else:
            photon.state = 0 
            
    elif chamber.ele[i].h > 0 and comp_ind == 2:
        if chamber.ele[i].end_type2 == 0:
            photon.state = 0 
            print('photon detected')
            ax.plot3D([Q[0], photon.O[0]], 
                      [Q[1], photon.O[1]], 
                      [Q[2], photon.O[2]], color = 'red')
        else:
            # photon.loc = i
            # for now assuming closed system, so detected or lamb reflects
            if np.random.rand() <= rhos:
                v, k = lamb_ref(photon.u, chamber.ele[i].n)
                photon.old.append(photon.O)
                ax.plot3D([Q[0], photon.O[0]], 
                          [Q[1], photon.O[1]], 
                          [Q[2], photon.O[2]], color = 'black')
                photon.O = Q
                photon.u = v
            else:
                photon.state = 0

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
                disp_r, *_ = disp_line_point(chamber.ele[i].P, chamber.ele[i].n, Q1)
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
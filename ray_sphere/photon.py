from distances import plane_create
import numpy as np

class photon_create:
    # convenient way to track the details of the photon
    def __init__(self, P, w, sigma_r, dt, state, loc): # o, u, v, state, loc, :  
        self.state = state
        # self.ref_status = ref_status
        self.loc = 0
        self.sigma_r = sigma_r
        self.dt = dt
        # self.P = P
        # self.w = w
        w1, w2 = plane_create(w)
        ### Need to choose if to you gaussian position spread or uniform
        self.O = P + w1*sigma_r*np.random.rand() + w2*sigma_r*np.random.rand()
        # self.o = c0 + w1*np.random.normal(0, sigma_r, 1)+ w2*np.random.normal(0, sigma_r, 1)
        self.u = w + w1*dt*np.random.rand() + w2*dt*np.random.rand()
        self.old = []

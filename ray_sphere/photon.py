from distances import plane_create
import numpy as np

class photon_create:
    info = []
    fate = []
    num = 0
    
    # convenient way to track the details of the photon
    def __init__(self, P, w, sigma_r, dt, state, loc=-99): # o, u, v, state, loc, :  
        photon_create.num += 1
        self.state = state
        self.loc = loc
        self.sigma_r = sigma_r
        self.dt = dt
        # self.P = P
        # self.w = w
        w1, w2 = plane_create(w)
        ### Need to choose if to you gaussian position spread or uniform
        self.O = P + w1*sigma_r*np.random.rand() + w2*sigma_r*np.random.rand()
        # self.o = c0 + w1*np.random.normal(0, sigma_r, 1)+ w2*np.random.normal(0, sigma_r, 1)
        self.u = w + w1*dt*np.random.rand() + w2*dt*np.random.rand()
    
    def __setattr__(self, name, value):
        self.__dict__[name] = value
        try:
            photon_create.info.append([photon_create.num, self.O, self.u, self.loc, self.state])
            # print(name, value)
        except:
            pass
        
    def info_rec(self, photon):
        photon_create.info.append([photon.num, photon.O, photon.u, photon.loc, photon.state])
    
    def info_empty(self):
        photon_create.info = []
    
    def fate_rec(self, photon, element, fate_type, angle=np.nan):
        # fate_type = {-1:"escaped", 0:"absorbed", 1:"detected"}
        photon_create.fate.append([photon.num, photon.O, photon.u, photon.loc, 
                                   element, fate_type, angle])
    
    def fate_empty(self):
            photon_create.fate = []
        
import numpy as np

def norm_e(v):
    """Euclidean norm of a."""
    return np.dot(v, v)**0.5

def Dist_sph(O, u, P0, r0):
    """Calculate the possible distances to the sphere for vector."""
    u = u/np.dot(u, u)**0.5 # Necessary?

    k = (np.dot(u, (O-P0))**2 - (np.dot(O-P0, O-P0) - (r0*r0)))**0.5
    if np.isnan(k):
        print("Doesn't intersect with sphere")
    else: # two options if the photon is in the sphere
        k1 = -np.dot(u, (O-P0)) - k # distance 1
        k2 = -np.dot(u, (O-P0)) + k # distance 2
        # Photon can move in +'ve and -'ve directions of u so assign 
        # v1 and a vector/pos. in direction of u.
        if np.dot(u, O + k1.u) >= 0: 
            dist1 = k1
            dist2 = k2
            v1 = O + k1*u
            v2 = O + k2*u
        else:
            dist2 = k1
            dist1 = k2
            v2 = O + k1*u
            v1 = O + k2*u
        return (v1, dist1, v2, dist2)

def line_line_int(P1, v1, P2, v2):
    # vector 1: P1 + t1.v1
    # vector 2: P2 + t2.v2
    # Position of intersection: Q = P1 + a.v1
    v1 = v1/norm_e(v1) # Necessary?
    v2 = v2/norm_e(v2) # Necessary?
    # cross of vt1, vt2, = 0, collinear
    v1xv2 = np.cross(v1 ,v2)
    if np.all(v1xv2==0):
        print('Collinear')
        if np.all(np.cross((P2-P1), v1)==0):
            print('Overlapping')
        else:
            print('Parallel')
    else: # lines will intersect
        # Find elements of v1xv2 that are nonzero and use the first
        # the ratio of an element of cross of P1, P2 with the correspondin 
        # element of v1xv2 gives the scalar that defines the intersection
        k = v1xv2 != 0 
        a = np.cross((P2-P1), v2)[k]/v1xv2[k]
        a = a[0]
        
        Q = P1 + a*v1
        return Q

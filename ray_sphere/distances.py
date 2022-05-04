import numpy as np

def norm_e(v):
    """Euclidean norm of a."""
    return np.dot(v, v)**0.5

def dist_sph(O, u, P0, r0):
    """Calculate the possible distances to the sphere for vector, returns
    Q1, dist1: intersection position and distance in direction of u
    Q2, dist2: other intersect, may also be in direction of u but greater
    """
    u = u/np.dot(u, u)**0.5 # Necessary?
    k = (np.dot(u, (O-P0))**2 - (np.dot(O-P0, O-P0) - (r0*r0)))**0.5
    if np.isnan(k):
        print("Doesn't intersect with sphere")
        Q1 = -1
        Q2 = -1
        dist1 = -1
        dist2 = -1
    else: # two options if the photon is in the sphere
        k1 = -np.dot(u, (O-P0)) - k # distance 1
        k2 = -np.dot(u, (O-P0)) + k # distance 2
        # Photon can move in +'ve and -'ve directions of u so assign 
        # v1 and a vector/pos. in direction of u.
        if np.dot(u, O + k1*u) >= 0: 
            dist1 = k1
            dist2 = k2
            Q1 = O + k1*u
            Q2 = O + k2*u
        else:
            dist2 = k1
            dist1 = k2
            Q2 = O + k1*u
            Q1 = O + k2*u
        return (Q1, dist1, Q2, dist2)

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
        return (Q, a)

def dist_line(O, P, p, u=[]): 
    # If u = [], dist is the shortest distance between O and P, 
    # If u is given, dist is the shortest distance between O and P along u

    if np.dot(O-P, O-P) == 0: 
        # O and P are at the same point, so distace = 0
        dist = 0
        D = O
    elif len(u) > 0:
        print('using line line interscept')
        (D, dist) = line_line_int(O, u, P, p)
        # Below is the old method but the line_line_int function can be used
    # elif len(u) > 0 and np.dot(p, u): # If u vector is given and moving towards line
    #     a = np.dot(P - O, p)/(np.dot(p, u))
    #     b = norm_e(P - O)
    #     dist = (a**2 - b**2)**0.5
    #     # Q = O + dist*u
    else:    
        # p = p/np.dot(p, p)**0.5 # ensure unit vector
        c = P + p*2*np.dot((O-P), (O-P))**0.5
        dist = np.dot(np.cross(c-P, P-O), np.cross(c-P, P-O))**0.5 \
            /np.dot(c-P, c-P)**0.5
        # Use pythagoras to calculate the point on the line 
        D = P + p*abs(np.dot(O-P, O-P) - dist**2)**0.5 \
            *np.dot(p, O-P)/(np.dot(p, O-P)*np.dot(p, O-P))**0.5
    return (dist, D)

def chord_dist(O, P, P0, r0):
    """
    Calculate the chord distance from photon on sphere to port.
    O: photon position
    P: port centre
    P0: sphere centre
    r0: sphere radius
    a: angle between the two points from the centre of the sphere
    dist: chord distance
    """

    a = np.arccos(np.dot(O-P0, P-P0) \
        /(np.dot(O-P0, O-P0)**0.5 * np.dot(P-P0, P-P0)**0.5))
    # chord = a*R
    dist = a*r0
    return dist

def plane_create(P0): 
    # Creates two orthogonal vectors, penpendicular to P0
    if (P0==0).all():
        print('No vector was provided')
    n = np.array([1, 0, 0])
    if (P0 == n).all():
        n = np.array([0, 1, 0])
    P0 = P0/np.dot(P0, P0)**0.5
    w2 = np.cross(P0,n)
    w2 = w2/np.dot(w2, w2)**0.5
    w1 = np.cross(w2,P0)
    if np.dot(w1,P0) < 0:
        w1 = np.cross(-w2,P0)
    w1 = w1/np.dot(w1, w1)**0.5
    return w1, w2  

def dist_cylinder(O, u, P, n, r): 
    # interactions with the cylinder walls not the caps
    # u is photon, n is port axis from c
    ## Update variables
    a = u - n*np.dot(u, n)
    b = (O-P) - n*np.dot(O-P, n)
    
    A = np.dot(a, a)
    B = 2*np.dot(a, b)
    C = np.dot(b, b) - r**2
    
    # Discriminant
    delta = B**2 - 4*A*C
    
    if delta < 0: 
        # doesn't interact with port, so create dummies
        s = -1
        v1 = -1
        d = -1
        v2 = -1
    else:            
        # 2 solutions of quadratic formula
        s1 = (-B + delta**0.5)/(2*A)
        s2 = (-B - delta**0.5)/(2*A)
        
        Q1 = O + s1*u
        # P2 = O + s2*u
        
        # from closest point on cylinder axis to intersection points dot u
        
        if np.dot((Q1 - (n*np.dot(Q1-P,n)+P)), u) > 0:
        # P2 - (n*np.dot(P2-c,n)+c)
            s = s1
            v1 = O + s*u
            d = s2 
            v2 = O + d*u
        else:
            s = s2
            v1 = O + s*u
            d = s1 
            v2 = O + d*u
    return (s, v1, d, v2)

# def dist_plane(O, u, P0, n):
#     # Distance to a plane with unit vector n that passes through pt c0
#     # if point c is on plane then n.(c - c0)
#     try: 
#         dist = (np.dot(n, P0) - np.dot(n, v1))/(np.dot(u, n))
#     except:
#         dist = -1
#     return dist
    

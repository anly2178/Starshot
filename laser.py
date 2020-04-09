"""
Here we define functions relating to the laser. These include:
    • The width of the laser
    • The average intensity
    • The fraction of power intercepted by the sail
"""

import numpy as np
from scipy import special

def find_width_at_L(L, wavelength, D):
    """
    Calculates the width of the laser beam, given the distance from the source,
    the diameter of the transmitter, and the laser wavelength.
    Using equation (38) from Kipping 2017.
    """
    W = 2**0.5 * L * wavelength / D
    return W

def find_avg_I_at_L(P, W):
    """
    Calculates the beam intensity at a distance L from the laser given the
    width W at a distance L from the laser array, and the power of the laser.
    Assumes laser dot is circular, beam width is rayleigh beam width,
    and intensity is uniform over circular dot.

    This gives the average intensity over the entire dot
    """
    I = 4*P/(np.pi * W**2)
    return I

def find_avg_I_over_sail(frac, P, A):
    """
    This gives average intensity over the sail.
    """
    I = frac*P / A
    return I

def find_frac_of_P_incident(A, W):
    """
    Calculates the fraction of laser power that is incident on the sail at a
    distance L from the laser array, given the area of the sail A and the width
    of the laser at distance L.
    Using equation (39) from Kipping 2017.
    """
    frac = (special.erf(np.sqrt(A/(2*np.pi))/W))**2
    return frac

P = 500 * 10**9
D = 1000
wavelength = 1200e-9
A = 10

L1 = 10000
L2 = 1000000000

W = find_width_at_L(L1,wavelength, D)
print(W)
frac = find_frac_of_P_incident(A, W)
print(frac)
print(find_avg_I_at_L(P, W))
print(find_avg_I_over_sail(frac, P, A))

#hello
W = find_width_at_L(L2,wavelength, D)
print(W)
frac = find_frac_of_P_incident(A, W)
print(frac)
print(find_avg_I_at_L(P, W))
print(find_avg_I_over_sail(frac, P, A))

import numpy as np
from scipy.special import erf

def find_beam_width(params, dist):
    """
    Finds beam width at some distance from the transmitter.
    Assumes circular Gaussian beam.
    """
    wavelength = params["wavelength"]
    diameter = params["diameter"]
    beam_width = np.sqrt(2) * dist * wavelength / diameter #m
    return beam_width

def find_fraction_incident(params, dist):
    """
    Finds fraction of laser power incident on the sail.
    Assumes circular Gaussian beam and projection of sail is circular.
        • Project for later: Find fraction incident on square sail.
    """
    #For small distances, produces zero division error because divides by small beam width.
    #Does not affect the calculations, so we ignore the error.
    with np.errstate(divide='ignore'):
        beam_width = find_beam_width(params, dist)
        area = params["area"]
        fraction = (erf(np.sqrt(area / (2*np.pi)) / beam_width))**2
    return fraction

# #This is from Kulkarni 2018; assumes beam is Bessel. Not relevant to Gaussian beam, and no longer compatible with library.
# def find_crit_dist(params):
#     """
#     Input:  The set of parameters
#     Output: The critical distance at which the laser begins to spill over the sail. (m)
#     """
#     #Get parameters
#     m_sail = params["m_sail"]
#     thickness = params["thickness"]
#     density = params["density"]
#     k = params["k"]
#     diameter = params["diameter"]
#     wavelength = params["wavelength"]
#     alpha = params["alpha"]
#
#     #Calculate the critical distance
#     sail_size = np.sqrt(m_sail / (k*density*thickness))
#     critical_dist = diameter * sail_size / (2*wavelength*alpha)
#     return critical_dist

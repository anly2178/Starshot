import numpy as np
from numpy import pi
from scipy.special import erf

def find_rayleigh_length(beam_waist, wavelength):
    """
    Returns rayleigh length (m) - distance from beam waist at which beam is equal to
    the sail size.
    """
    z_r = pi*beam_waist**2/wavelength #m
    return z_r


# def find_start_distance(params):
#     """
#     Returns the distance (m) away from the laser array that the sail should
#     be launched. At this distance the beam width is equal to the sail radius.
#     """
#     z_0 = find_focusing_length(params)
#     z_r = find_rayleigh_length(params)
#     start_d = z_0 - z_r
#     return start_d

# def find_acceleration_distance(params):
#     """
#     By assumption, we make acceleration distance the distance between the
#     two critical points where the size of the beam is equal to the size of
#     the sail. This is 2 times the rayleigh length.
#     """
#     return 2*find_rayleigh_length(params)

def find_crit_dist(params):
    """
    Assume dynamically varying focus of Gaussian beam. Beam is focused such that
    the beam waist is at the sail until a diffraction-limited distance, where
    the beam waist is equal to the sail radius.
    Returns the diffraction-limited distance (m).
    """
    D = 2*params["radius"] #(m) Diameter of sail
    d = params["diameter"] #(m) Diamter of transmitter array
    wavelength = params["wavelength"] #(m) Wavelength of laser
    L = pi*D*d/(4*wavelength)
    return L

def find_beam_width(params, dist):
    """
    Finds beam width at some distance from the transmitter.
    Assumes circular Gaussian beam, dynamically focused until diffraction-limited distance,
    where the beam spot is the size of the sail. Thereafter, beam waist is kept
    at that point.
    """
    crit_dist = find_crit_dist(params) #Distance where sail radius = beam width
    wavelength = params["wavelength"] #m
    diameter = params["diameter"] #m
    if dist <= crit_dist:
        beam_width = 2*wavelength*dist/(pi*diameter) #m
    elif dist > crit_dist:
        w_0 = 2*wavelength*crit_dist/(pi*diameter) #m
        d_waist = dist - crit_dist #m distance from beam waist
        z_r = find_rayleigh_length(w_0,wavelength) #m
        beam_width = w_0*(1+(d_waist/z_r)**2)**0.5 #m
    return beam_width

# def find_beam_width(params, dist):
#     """
#     Finds beam width at some distance from the transmitter.
#     Assumes circular Gaussian beam.
#     """
#     wavelength = params["wavelength"]
#     diameter = params["diameter"]
#     beam_width = np.sqrt(2) * dist * wavelength / diameter #m
#     return beam_width

def find_fraction_incident(params, dist):
    """
    Finds fraction of laser power incident on the sail.
    Assumes circular Gaussian beam and projection of sail is circular.
    """
    r = params["radius"]
    w = find_beam_width(params,dist)
    #To prevent overflow warnings, set a point where the fraction loss is appreciable; until then fraction = 1
    #Choose to care about fraction when it becomes <= 0.999, which occurs when:
    pt = 2*r**2/np.log(1000)
    if w**2 <= pt:
        fraction = 1
    else:
        fraction = 1-np.exp(-2*r**2/w**2)
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

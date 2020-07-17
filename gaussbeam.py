import numpy as np

def find_beam_width(diameter, wavelength, dist):
    """Calculate Gaussian beam width [m] at a distance from laser array."""
    beam_width = 2*wavelength*dist/(np.pi*diameter) #m
    return beam_width

def find_frac(radius, beam_width, dist):
    """Calculates fraction of laser power incidenet on circular lightsail."""
    #To prevent overflow warnings, set a point where the fraction loss is appreciable; until then fraction = 1
    #Choose to care about fraction when it becomes <= 0.999, which occurs when:
    pt = 2*radius**2/np.log(10000)
    if beam_width**2 <= pt:
        fraction = 1
    else:
        fraction = 1-np.exp(-2*radius**2/beam_width**2)
    return fraction

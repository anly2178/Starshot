import numpy as np
from scipy.special import erf

c = 2.998e8
stefan_boltzmann = 5.67037e-8

def find_beam_width(dist, wavelength, transmitter_D):
    """
    Finds beam width at some distance from the transmitter.
    Assumes circular Gaussian beam.
    """
    return np.sqrt(2) * dist * wavelength / transmitter_D #m

def find_fraction_incident(surface_area, beam_width):
    """
    Finds fraction of laser power incident on the sail.
    Assumes circular Gaussian beam and projection of sail is circular.
    """
    return (erf(np.sqrt(surface_area / 2*np.pi) / beam_width))**2

#Intensity is the derivative of power with respect to area - right?
#So if we find the fraction of power incident on a very small area, and
#divide by that very small area, we should have an approximation of intensity
#at that point. In particular, this will be the peak intensity.

def find_peak_intensity1(laser_power, beam_width):
    """
    This is using the equation for peak intensity. Can be found here:
    https://www.rp-photonics.com/gaussian_beams.html
    """
    peak_intensity = 2*laser_power / (np.pi * beam_width**2)
    return peak_intensity

def find_peak_intensity2(laser_power, surface_area, beam_width):
    """
    The fraction of laser power on a small area method.
    Should provide approximate peak intensity provided area is very small.
    """
    fraction = find_fraction_incident(surface_area, beam_width)
    incident_power = laser_power * fraction
    peak_intensity = incident_power/surface_area
    return peak_intensity

#Suppose:
dist = 1e7 #m
wavelength = 650e-9 #m
transmitter_D = 10 #m
laser_power = 500e9 #W
surface_area = 1e-11 #m^2

#Finding beam width
beam_width = find_beam_width(dist, wavelength, transmitter_D)

#Finding peak intensity using equation
peak_intensity1 = find_peak_intensity1(laser_power, beam_width)

#Finding peak intensity using the small area
peak_intensity2 = find_peak_intensity2(laser_power, surface_area, beam_width)

print("First method: {:e} Wm^-2 \nSecond method: {:e} Wm^-2".format(peak_intensity1,peak_intensity2))

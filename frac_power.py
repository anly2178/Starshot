import numpy as np
from scipy import special
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

np.seterr(divide='ignore', invalid='ignore')

c = 2.998e8
stefan_boltzmann = 5.67037e-8

def find_total_relative_energy(target_beta, reflectivity):
    """
    'Relative energy' is the energy of each photon expressed as a ratio
    of the sail's rest mass energy. Kipping 2017 eq 5
    Function calculates the relative energy summed over all incident photons.
    Assumes that photons are either absorbed or reflected, not transmitted.
    """
    beta = target_beta
    R = reflectivity
    return (-(1+R)*(1-beta)+np.sqrt(8*beta*R*(1-beta)+((1-beta)**2) * \
    ((1+R)**2)))/(4*R*(1-beta))

def find_firing_time(total_mass, surface_area, absorptivity, total_rel_energy,\
    max_temp, emissivity):
    """
    Finds the time taken for the sail to reach a target velocity,
    which is directly related to the total relative energy.
    Assumes the sail is in thermal equilibrium throughout firing time.
    """
    surface_density = total_mass / surface_area
    return (surface_density * absorptivity * (c**2) * total_rel_energy) / \
    (2 * emissivity * stefan_boltzmann * max_temp**4) #s

def find_acceleration_dist(target_beta, firing_time):
    """
    Finds distance travelled to reach 0.2c.
    Assumes constant acceleration.
    """
    return c * target_beta * firing_time / (1 + np.sqrt(1 - target_beta**2)) #m

def find_beam_width(dist, wavelength, transmitter_D):
    """
    Finds beam width at some distance from the transmitter.
    Assumes circular Gaussian beam.
    """
    return np.sqrt(2) * dist * wavelength / transmitter_D #m

def find_fraction_incident(surface_area, beam_width):
    """
    Finds fraction of laser power incident on the sail.
    Assumes circular Gaussian beam.
    """
    return (special.erf(np.sqrt(surface_area / 2*np.pi) / beam_width))**2

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

np.seterr(all='ignore', invalid='ignore')

c = 2.998e8
stefan_boltzmann = 5.67037e-8

def find_total_relative_energy(target_beta, reflectivity):
    """
    'Relative energy' is the energy of each photon expressed as a ratio
    of the sail's rest mass energy. Kipping 2017 eq 5
    Function calculates the relative energy summed over all incident photons.
    Assumes that photons are either absorbed or reflected, not transmitted.
    """
    B = target_beta
    R = reflectivity
    return (-(1+R)*(1-B)+np.sqrt(8*B*R*(1-B)+((1-B)**2) * ((1+R)**2)))/(4*R*(1-B))

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
    Assumes circular Gaussian beam and projection of sail is circular.
    """
    return (erf(np.sqrt(surface_area / 2*np.pi) / beam_width))**2

#
#The peak intensity is wrong for some reason. I essentially found a limit
#by considering a really small area. However, the peak intensity is such, that
#if you graphed the intensity distribution, the total power (area under the graph)
#is greater than the laser power.
# def find_peak_intensity(laser_power, dist, wavelength, transmitter_D):
#     """
#     Finds the peak intensity of the circular Gaussian beam.
#     """
#     beam_width = find_beam_width(dist, wavelength, transmitter_D)
#     peak_I = laser_power * find_fraction_incident(1e-10, beam_width)/1e-10
#     return peak_I
# #
# print(find_peak_intensity(500e9, 100000000, 650e-9,10))
#

#You could easily change this to show the intensity distribution if you have the peak intensity.
def plot_intensity(dist, wavelength, transmitter_D):
    """
    Plots the radial intensity/peak_intensity pattern at a distance away from the transmitter.
    Plots from x = -5 to x = 5 since the sail is likely not going to have a
    diameter of over 10m.
    Accepts distances at least 1000000 m, as otherwise it is not worthwhile considering
    diffraction effects.
    Function allows for user to overplot for multiple distances, however the user
    must call plt.show() after function call.
    """
    if dist < 1000000:
        print("Distance must be at least 1000 km for diffraction-limiting effects to be appreciable")
        return None

    beam_width = find_beam_width(dist, wavelength, transmitter_D)
    # peak_I = find_peak_intensity(laser_power, dist, wavelength, transmitter_D)

    radius = np.linspace(-5,5,1000)
    i = 0
    y = []
    while i < len(radius):
        r = radius[i]
        ratio_I = np.exp(-2 * r**2 / beam_width**2)
        y.append(ratio_I)
        i += 1
    plt.plot(radius,y)
    plt.title('Radial distribution of intensity')
    plt.ylabel(r'$\frac{I}{I_{peak}}$', fontsize = 16)
    plt.xlabel('Position along the diameter of beam (m)')
    return None

plot_intensity(1000000, 650e-9, 10)
plot_intensity(10000000, 650e-9, 10)
plot_intensity(30000000, 650e-9, 10)
plot_intensity(60000000, 650e-9, 10)
plot_intensity(100000000, 650e-9, 10)
plt.show()

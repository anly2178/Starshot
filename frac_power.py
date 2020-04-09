import numpy as np
from scipy import special
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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

"""
The following is for Kipping figure 6.
"""

mass_tot = 1e-3 #kg
sail_area = 16 #m^2 for one side
beta_targ = 0.2 #fraction of c
max_temp = 300 + 273.15 #K
wavelength = 650e-9 #m
emissivity = 1
transmitter_D = [1, 1e1, 1e2, 1e3] #m
absorptivity = np.logspace(-12,0,100) #absolute absorptivity
reflectivity = 1 - absorptivity #absolute reflectivity

fig = plt.figure()

i = 0
while i < len(transmitter_D):
    diameter = transmitter_D[i]
    r_tot = find_total_relative_energy(beta_targ, reflectivity)
    t = find_firing_time(mass_tot, sail_area, absorptivity, r_tot,\
    max_temp, emissivity)
    L = find_acceleration_dist(beta_targ, t)
    W = find_beam_width(L, wavelength, diameter)
    F = find_fraction_incident(sail_area, W)
    plt.plot(absorptivity, F)
    i += 1

plt.ylim(1e-10,1.5)
plt.ylabel('laser power fraction on the sail, $F_P$')
plt.yscale('log')

plt.xlim(1e-12, 1)
plt.xlabel('absorptivity of sail')
plt.xscale('log')


# def converse_functions(variables):
#     (A, t) = variables
#     beta = beta_targ
#     e = emissivity
#     surface_density = mass_tot / sail_area
#     eq = (surface_density * A * (c**2) * (-(1+1-A)*(1-beta)+\
#     np.sqrt(8*beta*(1-A)*(1-beta)+((1-beta)**2) * ((1+1-A)**2)))/\
#     (4*(1-A)*(1-beta))) / (2 * e * stefan_boltzmann * max_temp**4) - t
#     return [eq]
# surface_density = mass_tot / sail_area
# e = emissivity
# beta = beta_targ
# eq = lambda A : (surface_density * A * (c**2) * (-(1+1-A)*(1-beta)+\
# np.sqrt(8*beta*(1-A)*(1-beta)+((1-beta)**2) * ((1+1-A)**2)))/\
# (4*(1-A)*(1-beta))) / (2 * e * stefan_boltzmann * max_temp**4) - t
# A_guess = 1

# tick_times = np.array([1e-3, 1, 60, 3600, 3600 * 24, 3600 * 24 * 30])
# tick_location = np.array([])
# A_guesses = np.array([1e-10, 1e-7, 1e-6,1e-4,1e-3,1e-2])
#
# i = 0
# while i < len(tick_times):
#     t = tick_times[i]
#     eq = lambda A : (surface_density * A * (c**2) * (-(1+1-A)*(1-beta)+\
#     np.sqrt(8*beta*(1-A)*(1-beta)+((1-beta)**2) * ((1+1-A)**2)))/\
#     (4*(1-A)*(1-beta))) / (2 * e * stefan_boltzmann * max_temp**4) - t
#     A_guess = A_guesses[i]
#     solution = fsolve(eq, A_guess)
#     np.append(tick_location, solution)
#     i += 1

# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticklabels(["1 ms", "1 s", "1 min", "1 hr", "1 d", "1 mo"])
# ax2.set_xticks(tick_location)
# ax2.set_xlabel(r"continuous firing time")

plt.show()

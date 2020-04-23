import numpy as np
import time
from .TMM_analysis_sail.find_eq_temp import *
from .TMM_analysis_sail.optical_constants import n_silica


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
    r_tot = (-(1+R)*(1-B)+np.sqrt(8*B*R*(1-B)+((1-B)**2) * ((1+R)**2)))/(4*R*(1-B))
    return r_tot

def find_temp(params, beta, time):
    """
    Find the equilibrium temperature at a specific velocity and firing time.
    Using Kipping 2017 equation (37).
    """
    m_s = params["m_sail"]
    t = params["thickness"]
    rho = params["density"]
    R = params["reflectivity"]
    A = params["absorptance"]
    m_tot = 2*m_s
    area = m_s / (rho * t)
    sigma = m_tot / area #Effective surface density
    r_tot = find_total_relative_energy(beta, R)

    T = (r_tot * sigma * A * c**2 / (2 * stefan_boltzmann * time))**0.25
    return T

#The function below was experimental and now out of use.
"""
This is an adapted version of Justin's function, which originally solved
the equilibrium temperature for different absorption coefficients. It has been
adapted to solve equilibrium temperature for a given beta and absorption coefficient.
This function is SPECIFIC TO SILICON.
"""

def find_temp_silica(params, abs_coeff, beta):
    """
    This is an adapted version of Justin's function, which originally solved
    the equilibrium temperature for different absorption coefficients. It has been
    adapted to solve equilibrium temperature for a given beta and absorption coefficient.
    This function is SPECIFIC TO SILICA.
    """

    laser_power = params["power"] #W
    m_sail = params["m_sail"] * 1e3 #g
    ratio = laser_power/m_sail #W/g
    abs_coeff = [abs_coeff]

    thickness = params["thickness"] #m
    wavelength_0 = params["wavelength"] #m
    n = n_silica(wavelength_0)
    structure = [(n,-thickness)]

    density = params["density"]
    rho_S = density * thickness

    #Finding power absorbed, accouting for doppler shift
    wavelength = wavelength_0*np.sqrt((1+beta)/(1-beta))
    A = find_absorption_from_coefficient(structure, abs_coeff, wavelength)
    power_in = ratio*A*rho_S*(1-beta)/(1+beta)

    """ Note: Honestly, this below section should be made into its own function
        since it is a reusable block of code. Consider doing this at some point
        but for now, focus on optimising code and commenting
    """
    # The RHS is more complicated, since you can't get an expression for T explicitly
    # We need to integrate power flux over all wavelengths to get the total radiated power
    midpoint = 0
    bb_temp = (power_in/(2*1*5.67e-8))**0.25

    T_low = bb_temp             # Lower bound = max emissivity = black body temp
    T_high = bb_temp*10         # Upper bound arbitrary (might not hold at higher temps) - should find a way to set a true reasonable higher bound
    # Use trapezoidal rule to find the total power out for a given temperature
    def power_out(T):
        points = 101            # Can be changed if better resolution is required

        # Ilic paper uses 1-25 microns, but eqns should be valid from 0.5-50 microns if so required
        bounds = np.linspace(1e-6, 25e-6, points)
        power_out_at_wl = np.zeros(points)

        # Running each integral and adding to the list (optimisation here would be to fix list size and assign vals)
        i = 0
        for wavelength in bounds:
            power_out_at_wl[i] = (spectral_power_flux(wavelength, structure, T))
            i += 1
        power_out = np.trapz(power_out_at_wl, bounds, (25e-6-1e-6)/points)
        return power_out

    start_time = time.time()
    # Powers at the bounds of temperature interval
    P_high = power_out(T_high)
    P_low = power_out(T_low)

    # Halving the interval for a result
    while abs(P_high - P_low) >= 0.05*power_in:

    # The only issue we can really get is if P_high is too low - if this is
    # the case, just double P_high
        if (P_high <= power_in):
            T_high = T_high*2

        midpoint = (T_low+T_high)/2

        if power_out(midpoint) > power_in:
            T_high = midpoint
        else:
            T_low = midpoint

        P_high = power_out(T_high)
        P_low = power_out(T_low)
    # Take the midpoints as the final result since this is the result from halving interval
    midpoint = (T_high+T_low)/2
    print(midpoint)
    print("--- %s seconds ---" % (time.time() - start_time))
    return midpoint

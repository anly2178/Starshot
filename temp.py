import numpy as np
import time
from .TMM_analysis_sail.find_eq_temp import *
from .TMM_analysis_sail.optical_constants import n_silica
from .laser import find_fraction_incident

def find_max_power(params):
    """
    Returns a maximum power (W) that the laser can have for a given maximum temperature.
    """
    sb = 5.67e-8 #Stefan-Boltzmann constant
    T_max = params["max_temp"] #K, max temperature sail can sustain
    A = params["area"] #m^2, area of sail
    absorb = params["absorptance"] #Absolute absorption of sail
    #Assume max temperature occurs at beta=0, dist=0
    beta = 0
    dist = 0
    #Bounds for halving the interval
    P_bb = 2*A*sb*T_max**4 #W, Power incident on black body
    P_high = P_bb/absorb #W, Upper bound; power incident on sail accounting for absorptance
    P_low = P_high/1000 #W, Arbitrary lower bound
    #Create structures to not accidentally destroy original parameters
    params_high = params.copy()
    params_high["power"] = P_high
    params_low = params.copy()
    params_low["power"] = P_low
    #Temperatures at bounds
    T_high = find_one_temp(params_high,beta,dist)
    T_low = find_one_temp(params_low,beta,dist)
    #Midpoint temperature
    T_mid = (T_high+T_low)/2
    #Check the run time
    start_time = time.time()
    while abs(T_mid - T_max) >= 0.001*T_max:
        #If the high temperature is too low for some reason
        if T_high <= T_max:
            params_high["power"] = 2*params_high["power"]
        #If the low temperature is too high for some reason
        if T_low >= T_max:
            params_low["power"] = params_low["power"]/2
        #Midpoint
        P_high = params_high["power"]
        P_low = params_low["power"]
        P_mid = (P_high + P_low)/2
        params_mid = params.copy()
        params_mid["power"] = P_mid
        T_mid = find_one_temp(params_mid,beta,dist)
        if T_mid > T_max:
            params_high["power"] = P_mid
        else:
            params_low["power"] = P_mid
        print(P_mid)
    print("--- %s seconds ---" % (time.time() - start_time))
    P_mid = (params_high["power"] + params_low["power"])/2
    return P_mid

def find_one_temp(params, beta, dist):
    """
    Returns the equilibrium temperature at a specific speed and distance,
    accounting for doppler shift and diffraction effects.

    Input:  parameters, beta (fraction of c), distance from DE system (m)
    Ouput:  Equilibrium temperature (K)
    """

    #Extracting the parameters
    material = params["material"]
    laser_power = params["power"] #W
    m_sail = params["m_sail"] * 1e3 #g
    ratio = laser_power/m_sail #W/g
    abs_coeff = [params["abs_coeff"]]
    thickness = params["thickness"] #m
    wavelength_0 = params["wavelength"] #m
    area = params["area"] #m^2
    rho_S = m_sail / area #gm^-2

    #Determining which material to calculate absorption for
    if material == 'SiO2':
        n = 1.45
    elif material == 'GeO2':
        n = 1.6
    #Create structure
    i = 0
    structure = []
    while i < len(thickness):
        if i%2 == 0:
            layer = (n,-thickness[i]) #Material
        elif i%2 == 1:
            layer = (1,-thickness[i]) #Vacuum
        structure.append(layer)
        i += 1
    #Finding power absorbed, accounting for doppler shift and diffraction effects
    wavelength = wavelength_0*np.sqrt((1+beta)/(1-beta))
    A = find_absorption_from_coefficient(structure, abs_coeff, wavelength)
    try:
        fraction = find_fraction_incident(params, dist)
    except TypeError: #Diameter not given, which occurs when finding max power
        fraction = 1 #Very good approx for small distances
    power_in = ratio*A*rho_S*(1-beta)/(1+beta) * fraction

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
            poawl = (spectral_power_flux(wavelength, structure, T))
            if poawl == None: #In the case of overflow error, return None
                return None
            else:
                power_out_at_wl[i] = poawl
            i += 1
        power_out = np.trapz(power_out_at_wl, bounds, (25e-6-1e-6)/points)
        return power_out

    start_time = time.time()
    # Powers at the bounds of temperature interval
    P_high = power_out(T_high)
    P_low = power_out(T_low)
    if P_high == None or P_low == None: #If overflow error, return temperature as 0
        return 0

    # Halving the interval for a result
    while abs(P_high - P_low) >= 0.01*power_in:
    # The only issue we can really get is if P_high is too low - if this is
    # the case, just double P_high
        if (P_high <= power_in):
            T_high = T_high*2

        midpoint = (T_low+T_high)/2
        P_mid = power_out(midpoint)
        if P_mid == None: #If overflow error, return temperature as 0
            return 0
        if P_mid > power_in:
            T_high = midpoint
        else:
            T_low = midpoint

        P_high = power_out(T_high)
        P_low = power_out(T_low)
        if P_high == None or P_low == None: #If overflow error, return temperature as 0
            return 0
        print(midpoint)
    # Take the midpoints as the final result since this is the result from halving interval
    midpoint = (T_high+T_low)/2
    print("--- %s seconds ---" % (time.time() - start_time))
    return midpoint


def add_all_temp(params, state):
    """
    Appends an array of equilibrium temperatures, corresponding to the state of
    the sail at each time (state includes speed and distance).

    Input:  parameters (dictionary)
            state (multidimensional array)
    Output: The multidimensional array of states with the equilibrium temperatures
            appended as the third column.

    Example input:      state = [[beta1,beta2,...,betan],
                                  [dist1,dist2,...,distn]]

    Example output:     state = [[beta1,beta2,...,betan],
                                  [dist1,dist2,...,distn],
                                  [temp1,temp2,...,tempn]]
    """

    #Extract speeds and distances
    betas = state[0,:]
    dists = state[1,:]
    n = len(betas)

    #Create structure holding temperatures
    temps = np.zeros(n)

    i = 0
    while i < n:
        beta = betas[i]
        dist = dists[i]
        temp = find_one_temp(params, beta, dist)
        temps[i] = temp
        i += 1
    new_state = np.vstack((state, temps))
    return new_state

#The below functions underestimate the equilibrium temperature as they assume
#the sail is a blackbody. To provide more accurate results, we opt for
#functions that are material/structure specific such as the ones above.

# c = 2.998e8
# stefan_boltzmann = 5.67037e-8
#
# def find_total_relative_energy(target_beta, reflectance):
#     """
#     'Relative energy' is the energy of each photon expressed as a ratio
#     of the sail's rest mass energy. Kipping 2017 eq 5
#     Function calculates the relative energy summed over all incident photons.
#     Assumes that photons are either absorbed or reflected, not transmitted.
#     """
#     B = target_beta
#     R = reflectance
#     r_tot = (-(1+R)*(1-B)+np.sqrt(8*B*R*(1-B)+((1-B)**2) * ((1+R)**2)))/(4*R*(1-B))
#     return r_tot

# def find_temp(params, beta, time):
#     """
#     Find the equilibrium temperature at a specific velocity and firing time.
#     Using Kipping 2017 equation (37).
#     """
#     m_s = params["m_sail"]
#     t = params["thickness"]
#     rho = params["density"]
#     R = params["reflectance"]
#     A = params["absorptance"]
#     m_tot = 2*m_s
#     area = m_s / (rho * t)
#     sigma = m_tot / area #Effective surface density
#     r_tot = find_total_relative_energy(beta, R)
#
#     T = (r_tot * sigma * A * c**2 / (2 * stefan_boltzmann * time))**0.25
#     return T

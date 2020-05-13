import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from sympy.solvers import solve
from sympy import Symbol
from .tmm import tmm, general_tmm
from .make_transfer_matrix import make_transfer_matrix
from .optical_constants import n_silica, n_germania
import time
import warnings

" ============================================================================ "

""" Finds the absolute (hemispherical) absorption from the coefficient given
    This function should be able to apply to any structure given it only has
    a real portion of RI in structure.

    structure: [(n,d), (n,d)], each tuple relates to own layer. d in m
    abs_coeff: Given as a list/tuple, each value in cm-1
    wavelength: Incident wavelength (of laser), given in m

    Output: List of absolute absorption, calculated using TMM at *NORMAL INCIDENCE*
"""

def find_absorption_from_coefficient(structure, abs_coeff, wavelength):
    # This block of code finds the absorption given the absorption coefficient
    R = np.array([None]*len(abs_coeff))
    T = np.array([None]*len(abs_coeff))

    k = 0
    while k < len(abs_coeff):
        im_RI = 1j*wavelength*100*abs_coeff[k]/(4*pi)
        # need to create temp_struc to be the same size as structure so it can be filled
        # without accidentally overwriting structure, so we do this:
        temp_struc = [None]*len(structure)
        for l in range(len(structure)):
            temp_struc[l] = [None]*2    # since each tuple is of form (n,d) only
        i = 0
        while i < len(structure):
            if structure[i][0] != 1:
                temp_struc[i][0] = structure[i][0] + im_RI   # structure in this iteration
            else:
                temp_struc[i][0] = structure[i][0]
            temp_struc[i][1] = structure[i][1]      # width isnt going to change
            i += 1

        optical_parameters = tmm(temp_struc, wavelength)
        R[k] = (optical_parameters[0]*np.conj(optical_parameters[0])).real
        T[k] = (optical_parameters[1]*np.conj(optical_parameters[1])).real
        k += 1

    A = [1]*(len(abs_coeff)) - R - T
    return A

" ============================================================================ "

""" This function simply defines (hemispherical) emissivity as a function of
    theta so it can be integrated (directional emissivity). Use Kirchoff's
    Law - only need to find absorption as a function of theta (incident angle).
    Note this assumes emissivity does not vary around the surface of the sail
    (perfectly smooth?). Absorption = 1 - reflectance - transmittance.

    Structure required for thickness parameters, wavelength used to calculate
    optical constants (using equations set by papers)

    theta: incident angle in radians
    structure: [(n,d), (n,d)], each tuple relates to own layer. d in m
    wavelength: given in m

    Output: *NOT ACTUALLY DIRECTIONAL_EMISSIVITY* - The output includes a
            cos(theta)*sin(theta) term multiplied to the directional emissivity
            which is present in the Jacobian when changing variables for an
            integral to find the spectral hemispherical emissivity.
            Otherwise, it is basically the directional emissivity evaluated
            at a specific wavelength
"""

    # NOTE: WE ARE CURRENTLY LOOKING AT SILICON AND SILICA MOSTLY, SO IT IS
    # OKAY TO USE OPTICAL CONSTANTS FROM STRUCTURE TO IDENTIFY THEM TO FIND
    # OPTICAL CONSTANTS

def directional_emissivity(theta, structure, wavelength):
    # First, we need to find the optical constants using wavelength

    n1 = n_silica(wavelength)
    n2 = n_germania(wavelength)

    # Section of code that creates a new structure list that is based on
    # calculated optical constants using wavelength
    temp_struc = [None]*len(structure)
    for l in range(len(structure)):
        temp_struc[l] = [None]*2    # since each tuple is of form (n,d) only
        # Assigning optical constants as we go
        if structure[l][0] == 1.45:
            temp_struc[l][0] = n1
        elif structure[l][0] == 1.6:
            temp_struc[l][0] = n2
        else:
            temp_struc[l][0] = structure[l][0]
        # Keeping the thickness of each layer
        temp_struc[l][1] = structure[l][1]

    # This block gives an expression for emissivity in terms of theta (and wavelength)
    # First set out by finding the reflectance and transmittance of structure at
    # this wavelength and angle
    optical_parameters = general_tmm(temp_struc, wavelength, theta)
    R = (optical_parameters[0]*np.conj(optical_parameters[0]) + optical_parameters[2]*np.conj(optical_parameters[2]))/2
    T = (optical_parameters[1]*np.conj(optical_parameters[1]) + optical_parameters[3]*np.conj(optical_parameters[3]))/2
    dEpsilon = (1-R-T)*cos(theta)*sin(theta)
    return dEpsilon.real

" ============================================================================ "

""" Finds the spectral power flux of a pseudo-real (perfectly flat and smooth)
    sail that is SYMMETRICAL on both sides of the sail. Look towards making
    changes to this function or perhaps a new function to account for asymmetry.

    Considers inputs to return power flux at specific wavelengths and temperatures

    wavelength: given in m
    structure: [(n,d), (n,d)], each tuple relates to own layer. d in m
    temperature: given in K
"""

    # Here, we need an expression in terms of the wavelength since we need to
    # integrate the spectral emissivity. An integral should also be here to convert
    # the directional_emissivity to hemispherical emissivity

def spectral_power_flux(wavelength, structure, temperature):
    # First, give expression for radiation emitted by a black body
    h = 6.62607004e-34       # Planck's constant in SI
    c = 299792458             # speed of light in SI
    k_B = 1.38064852e-23        # Boltzmann constant in SI
    sigma = 5.67e-8
    with warnings.catch_warnings(): #Overflow error occurs at low temperatures 
        warnings.filterwarnings('error')
        try:
            I = ((2*h*c**2)/wavelength.astype(float)**5)*(1/(np.exp(h*c/(wavelength.astype(float)*k_B*temperature.astype(float)))-1))         # Planck's Law
        except Warning as e:
            print('Error found:', e)
            return None

    # Now give expression for hemispherical emissivity. Note factor of 4: 2 comes
    # from integrating wrt phi, the other to account for both faces of sail

    # Use trapezoidal integration to speed things up
    points = 50         # Number of points in the integration
    bounds = np.linspace(0,pi/2,points)
    spec_hemi_ems = points*[None]
    i = 0
    for theta in bounds:
        spec_hemi_ems[i] = (4*directional_emissivity(theta, structure, wavelength))
        i += 1
    power_flux = pi*I*np.trapz(spec_hemi_ems, bounds, pi/2/points)
    return power_flux


" ============================================================================ "

""" Calculate the equilibrium temperature at given absorption coefficients.
    Requires following args: ratio (P/m ratio) in W/g, structure, rho_S (g/m^2)

    Will return three lists which give the midpoint and interval values for
    equilibrium temperature solved by halving the interval at each
    absorption coefficient
"""

def solveTemp(ratio, structure, rho_S, abs_coeff = np.logspace(-2,-6,15)):
    # Absorbance is measured as 1-R-T where R, T are reflectance and transmittance.
    # dependence on the absorption coefficient is based on definition of RI from paper.

    """ TURNS OUT, WE HAVE TO ACCOUNT FOR THE DOPPLER SHIFT!
    THE ILIC PAPER PLOTS THE MAXIMUM TEMPERATURE CURVE BASED ON THE SPEED OF
    THE CRAFT. HENCE WE HAVE TO FIND THE VALUE FOR BETA WHERE THE POWER IN PER
    UNIT AREA IS MAXIMISED AS WELL. SO WE RUN THIS THROUGH CODE BELOW:
    """
    wavelength_0 = 1.2e-6
    beta = np.linspace(0,0.2,100)
    power_in = [0]*10       # need to find maximum p_in based on beta

    # Loop that gets the maximum power value (change to vector to optimise?)
    for b in beta:
        wavelength = wavelength_0*np.sqrt((1+b)/(1-b))
        A = find_absorption_from_coefficient(structure, abs_coeff, wavelength)
        # Finding the LHS of Ilic equation
        power_beta = ratio*A*rho_S*(1-b)/(1+b)

        if power_beta[-1] > power_in[-1]:
            power_in = power_beta

    """ Note: Honestly, this below section should be made into its own function
        since it is a reusable block of code. Consider doing this at some point
        but for now, focus on optimising code and commenting
    """
    # The RHS is more complicated, since you can't get an expression for T explicitly
    # We need to integrate power flux over all wavelengths to get the total radiated power
    temps = []
    midpoints = []
    highs = []
    lows = []
    for P in power_in:              # Related to each x-val (abs coeff)
        start_time = time.time()
        bb_temp = (P/(2*1*5.67e-8))**0.25

        T_low = bb_temp             # Lower bound = max emissivity = black body temp
        T_high = bb_temp*10         # Upper bound arbitrary (might not hold at higher temps) - should find a way to set a true reasonable higher bound

        # Use trapezoidal rule to find the total power out for a given temperature
        def power_out(T):
            points = 101            # Can be changed if better resolution is required

            # Ilic paper uses 1-25 microns, but eqns should be valid from 0.5-50 microns if so required
            bounds = np.linspace(1e-6, 25e-6, points)
            power_out_at_wl = points*[None]

            # Running each integral and adding to the list (optimisation here would be to fix list size and assign vals)
            i = 0
            for wavelength in bounds:
                power_out_at_wl[i] = (spectral_power_flux(wavelength, structure, T))
                i += 1
            power_out = np.trapz(power_out_at_wl, bounds, (25e-6-1e-6)/points)
            return power_out

        # Powers at the bounds of temperature interval
        P_high = power_out(T_high)
        P_low = power_out(T_low)

        # Halving the interval for a result
        while abs(P_high - P_low) >= 0.05*P:

        # The only issue we can really get is if P_high is too low - if this is
        # the case, just double P_high
            if (P_high <= P):
                T_high = T_high*2

            midpoint = (T_low+T_high)/2

            if power_out(midpoint) > P:
                T_high = midpoint
            else:
                T_low = midpoint

            P_high = power_out(T_high)
            P_low = power_out(T_low)

        # Take the midpoints as the final result since this is the result from halving interval
        midpoints.append((T_high+T_low)/2)

        # Also keep interval bounds in case need to compare (maybe also used to check error/give interval)
        highs.append(T_high)
        lows.append(T_low)

        # Timer and printing out the midpoint temperature in case needs to be seen
        print(T_high/2+T_low/2)
        print("--- %s seconds ---" % (time.time() - start_time))

    temps = [midpoints, highs, lows]

    return temps


" ============================================================================ "


























# Useless for now


def solveTemp_csi(ratio, structure, rho_S, abs_coeff = np.logspace(1e-2,10^-13,1000), wavelength = 1.2e-6, eps = 0.4):
    # Absorbance is measured as 1-R-T where R, T are reflectance and transmittance.
    # dependence on the absorption coefficient is based on definition of RI from paper.

    R = np.array([None]*1000)
    T = np.array([None]*1000)

    k = 0
    while k < len(abs_coeff):
        im_RI = 1j*wavelength*100*abs_coeff[k]/(4*pi)
        # need to create temp_struc to be the same size as structure so it can be filled
        # without accidentally overwriting structure
        temp_struc = [None]*len(structure)
        for l in range(len(structure)):
            temp_struc[l] = [None]*2    # since each tuple is of form (n,d) only
        i = 0
        while i < len(structure):
            if structure[i][0] == 3.5:
                temp_struc[i][0] = structure[i][0] + im_RI   # structure in this iteration
            else:
                temp_struc[i][0] = structure[i][0] + 1j*wavelength*100*1e-6/(4*pi)
            temp_struc[i][1] = structure[i][1]      # width isnt going to change
            i += 1
        R[k] = (tmm(temp_struc, wavelength)[0]*np.conj(tmm(temp_struc, wavelength)[0])).real
        T[k] = (tmm(temp_struc, wavelength)[1]*np.conj(tmm(temp_struc, wavelength)[1])).real
        k += 1

    A = [1]*1000 - R - T

    power_in = ratio*A*rho_S

    sigma = 5.67e-8
    T0 = 2.73
    T = (power_in/(2*eps*sigma)+T0**4)**0.25

    return T

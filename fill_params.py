import numpy as np
from numpy import sin, cos, pi
import scipy
import scipy.integrate as integrate
from .TMM_analysis_sail.optical_constants import n_silica, n_germania
from .TMM_analysis_sail.tmm import tmm

"""
Fills in missing parameters from the given ones, if possible.
"""

def fill_thickness(params):
    """
    Input:  Set of parameters that contains
                • mass of sail
                • density
                • surface area
    Output: Dictionary of parameters with thickness (m) added
    """
    m_sail = params["m_sail"] #kg
    density = params["density"] #kg/m^3
    area = params["area"]

    thickness = m_sail / (density * area)
    params["thickness"] = thickness
    return params

def fill_area(params):
    """
    Input:  Set of parameters that contains
                • mass of sail
                • density
                • thickness
    Output: Dictionary of parameters with surface area (m^2) added
    """
    m_sail = params["m_sail"] #kg
    density = params["density"] #kgm^-3
    thickness = params["thickness"] #m

    area = m_sail / (density * thickness) #m^2
    params["area"] = area
    return params

def fill_radius(params):
    """
    Input:  Set of parameters that contains
                • mass of sail
                • density
                • thickness
    Output: Dictionary of parameters with radius of circular sail (m) added
    """
    area = params["area"]
    radius = (area / pi)**0.5
    params["radius"] = radius
    return params

def find_structure(params):
    """
    Finds the structure - for the next two functions.
    """
    thickness = params["thickness"]
    material = params["material"]
    wavelength = params["wavelength"]
    #Find the refractive index of the material at the wavelength
    if material == 'SiO2':
        n = n_silica(wavelength)
    elif material == 'GeO2':
        n = n_germania(wavelength)
    else:
        print("""List of materials available and their corresponding strings:
                Material              String
                Germania ------------ 'GeO2'
                Silica   ------------ 'SiO2'""")
        return None

    #Set up arguments for tmm
    structure = [(n,-thickness)]
    return structure


def fill_abs_ref_tra(params):
    """
    Input:  Dictionary of parameters, that must include the matrial and absorption coefficient
            and must have 'None' entered for absorptance, reflectance and transmittance
    =====================================================================
    List of materials and their corresponding strings
        Material              String
        Germania ------------ 'GeO2'
        Silica   ------------ 'SiO2'
        Will be updated
    =====================================================================
    Output: Dictionary of parameters with absorptance, reflectance and transmittance
            included.
    """

    #Extract parameters
    wavelength = params["wavelength"]
    abs_coeff = params["abs_coeff"]
    #Finds the structure
    structure = find_structure(params)
    if structure == None:
        return params
    abs_coeff = [abs_coeff]

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

    params["absorptance"] = A[0]
    params["reflectance"] = R[0]
    params["transmittance"] = T[0]
    return params

def fill_W(params):
    """
    Fills the square root of RAAD as defined in Ilic 2018.
    Units: sqrt(g)/m
    """
    #Set up problem
    if not "angles_coeffs" in params:
        structure = find_structure(params)
        if structure == None:
            return params
    elif "angles_coeffs" in params:
        structure = params["reflectance"]
    area = params["area"] #m^2
    m_sail = params["m_sail"]*1000 #g
    wavelength = params["wavelength"] #m
    rho_S = m_sail / area #Surface density g/m^2
    #To evaluate the integral we define dW
    def dW(beta, structure, rho_S, wavelength):
        gamma = 1/np.sqrt(1-beta**2)
        ds_wavelength = wavelength*np.sqrt((1+beta)/(1-beta))
        if not "angles_coeffs" in params:
            dW = np.sqrt(rho_S)/(tmm(structure, ds_wavelength)[0]*np.conj(tmm(structure, ds_wavelength)[0])).real * (gamma*beta)/(1-beta)**2
        elif "angles_coeffs" in params:
            dW = np.sqrt(rho_S)/structure * (gamma*beta)/(1-beta)**2
        return dW

    target = params["target"]
    W = integrate.quad(dW, 0, target, args=(structure, rho_S, wavelength))
    params["W"] = W[0]
    return params

def fill_power(params):
    """
    Fills the power needed to reach the target speed at a target acceleration
    distance.
    Params must contain a target acceleration distance.
    """
    c = 2.998e8 #m/s
    m = params["m_sail"]*1000 #g
    W = params["W"] #sqrt(g)/m
    D = params["radius"]*2 #m, diameter of sail
    L = params["accel_dist"] #m,
    P = c**3*np.sqrt(m)*W*D/(1000*L) #W
    params["power"] = P
    return params

def fill_diameter(params):
    """
    If None is given as diameter, fill the diameter with the same diameter
    as the sail. This means the sail will start from the laser array - for Gaussian beam.
    For Ilic beam, diameter is constrained by other parameters.
    (m)
    """
    c = 2.998e8
    wavelength = params["wavelength"]
    m_sail = params["m_sail"] * 1000 #g
    W = params["W"]
    P = params["power"] #GW
    d = (2*wavelength*c**3*np.sqrt(m_sail)*W)/(1000*P)
    params["diameter"] = d
    return params

def fill_params(params):
    """
    Fills the missing parameters if possible.
    """
    for key, value in params.items():
        if value == None:
            if key == 'thickness':
                fill_thickness(params)
            elif key == 'area':
                fill_area(params)
            elif key == 'radius':
                fill_radius(params)
            elif key == 'absorptance' or key == 'reflectance' or key == 'transmittance':
                fill_abs_ref_tra(params)
            elif key == 'W':
                fill_W(params)
            elif key == 'diameter':
                fill_diameter(params)
    return params

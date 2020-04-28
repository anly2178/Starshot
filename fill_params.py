import numpy as np
from numpy import sin, cos, pi
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
    m_sail = params["m_sail"]
    density = params["density"]
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
    m_sail = params["m_sail"]
    density = params["density"]
    thickness = params["thickness"]

    area = m_sail / (density * thickness)
    params["area"] = area
    return params

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
    thickness = params["thickness"]
    abs_coeff = params["abs_coeff"]
    material = params["material"]

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
            elif key == 'absorptance' or key == 'reflectance' or key == 'transmittance':
                fill_abs_ref_tra(params)
    return params

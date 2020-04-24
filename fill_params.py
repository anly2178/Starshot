from .TMM_analysis_sail.find_eq_temp import *
from .TMM_analysis_sail.optical_constants import n_silica, n_germania

"""
Fills in missing parameters from the given ones, if possible.
"""

def add_thickness(params):
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

def add_area(params):
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

"""
The following function is used if one has:
    1. ABSORPTION COEFFICIENT of
    2. MATERIAL that is either:
            • silica
            • germania
            • (This will continue to be updated)

If material is not available in the list, one must include their own value
of absorptance into their set of parameters.
"""

def add_absorptance(params):
    """
    Input:  material is a STRING representing the material of the sail.
            thickness of the sail (m)
            absorption coefficient (cm^-1)
            wavelength of laser (m)
    =====================================================================
    List of materials and their corresponding strings
        Material              String
        Germania ------------ 'GeO2'
        Silica   ------------ 'SiO2'
        Will be updated
    =====================================================================
    Output: Absorptance - fraction of intensity absorbed. (absolute absorption)

    The list of parameters should be defined AFTER absorptance is found.
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

    structure = [(n,-thickness)]
    abs_coeff = [abs_coeff]

    absorptance = find_absorption_from_coefficient(structure, abs_coeff, wavelength)
    absorptance = absorptance[0]
    params["absorptance"] = absorptance
    return params

def fill_params(params):
    """
    Fills the missing parameters if possible.
    """
    for key, value in params.items():
        if value == None:
            if key == 'thickness':
                add_thickness(params)
            elif key == 'area':
                add_area(params)
            elif key == 'absorptance':
                add_absorptance(params)
    return params

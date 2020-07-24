import os
import pickle
from save_load_mat import has_saved
from material_helpers import interpolate_from_list, use_equation

""" Each material should have a:
        - name (str; potentially multiple names. Main name is chemical formula
          e.g. SiO2, GeO2, Si3N4, etc.)
        - density (float; kg/m^3)
        - set of optical constants
            * real component (list of list/tuples, each tuple has 2 parts,
              the wavelength and the value)
            * imaginary component (same as above)
        - boolean value to designate whether it has equations which can describe
          its optical constants. If this is the case, the get_n or get_k
          functions will be altered in order to use the equations. For an equation
          to be used, the user must code it in manually
            (see materials_equations.py)
"""

class Material:

    def __init__(self, name, density, n_list = None, k_list = None, has_equations_for_n = False, has_equations_for_k = False):
        """ Constructor requires at least the name and the density
        """
        self.name = name
        self.density = density
        self.n_list = n_list
        self.k_list = k_list
        self.has_equations_for_n = has_equations_for_n
        self.has_equations_for_k = has_equations_for_k
        if not has_saved(self):
            self._save_material()

    def get_density(self):
        return self.density

    def set_density(self, density):
        self.density = density

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def make_list_from_file(path):
        """ Takes in the file path of a CSV with each entry organised as
            [wavelength],[n/k] and forms a list of tuples. If any tuples are a
            different size to any of the others, will raise a ValueError.
            If the file is not a CSV, will also raise a ValueError.

            Perhaps future implementation for tab separated files?
        """
        f = open(absolute_path, 'r')
        entries = f.readlines()
        f.close()
        for entry in entries:
            float(entry.strip().split(','))
        list.sort(entries)      # Sort in ascending order in terms of wavelength
        return entries

    def set_n_list(self, absolute_path):
        """ Takes in a string which is the ABSOLUTE (relative should work too,
            but just to be safe) path of the file that contains the data for n.

            File should be a CSV with each entry organised as:
            [wavelength],[n]

            If any of the entries have a length that is not 2, raises a ValueError.
            WAVELENGTHS ARE IN MICROMETRES (this is the most convenient here,
            and it seems refractiveindex.info tends to record in microns as well
            if any data is taken from there)
        """
        self.n_list = make_list_from_file(absolute_path)

    def get_n(self, wavelength):
        """ If an equation needs to be used, it will use an equation. Each
            equation in materials_equations is identified by the material name
            (and if there is an equation for n and k, the 'n' and 'k' string
            identifiers in the second argument helps to differentiate).
            Requires wavelength as float to find values
        """
        if has_equations_for_n:
            n = use_equation(self, 'n', wavelength)
        else:
            n = interpolate_from_list(self.n_list, wavelength)
        return n

    def get_n_list(self):
        return self.n

    def set_k_list(self, absolute_path):
        """ Same as above, but for k
        """
        self.k_list = make_list_from_file(absolute_path)

    def get_k(self, wavelength):
        if has_equations_for_k:
            k = use_equation(self, 'k', wavelength)
        else:
            k = interpolate_from_list(self.k_list, wavelength)
        return k

    def get_k_list(self):
        return self.k_list

    def get_has_equations_for_n(self):
        return has_equations_for_n

    def set_has_equations_for_n(self, bool):
        has_equations_for_n = bool

    def get_has_equations_for_k(self):
        return has_equations_for_n

    def set_has_equations_for_k(self, bool):
        has_equations_for_n = bool

    def _save_material(self):
        with open('material_data.pkl', 'ab') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    # def save_material(self):
    #    """Saves the material properties to a file for ease of use in the future"""
    #    TO DO
    #
    # def load_material(material_name):
    #    """Loads the material properties from a file based on the material_name (string) given"""
    #    TO DO




    # Possible implementation in future?
    # def set_optical_constants(self, absolute_path):
    #     """ Takes in a string which is the ABSOLUTE path of the file that contains
    #         the data for BOTH n AND k. If the user has a separate file for each,
    #         please see methods:
    #             - set_n(self, absolute_path)
    #             - set_k(self, absolute_path)
    #
    #         This file should be a CSV so that each entry is organised as:
    #         [wavelength],[n],[k]
    #
    #         If any of the entries have a length that is not 3, raises a
    #         ValueError. Ensure the wavelength is in
    #     """

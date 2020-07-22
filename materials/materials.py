import os
from materials_equations import use_equation

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

    def __init__(self, name, density, n = None, k = None has_equations = False):
        """ Constructor requires at least the name and the density
        """
        self.name = name
        self.density = density
        self.n = n
        self.k = k
        self.has_equations = has_equations

    def get_density(self):
        return self.density

    def set_density(self, density):
        self.density = density

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

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
        return entries

    def set_n(self, absolute_path):
        """ Takes in a string which is the ABSOLUTE (relative should work too,
            but just to be safe) path of the file that contains the data for n.

            File should be a CSV with each entry organised as:
            [wavelength],[n]

            If any of the entries have a length that is not 2, raises a ValueError.
            WAVELENGTHS ARE IN MICROMETRES (this is the most convenient here,
            and it seems refractiveindex.info tends to record in microns as well
            if any data is taken from there)
        """
        self.n = make_list_from_file(absolute_path)

    def get_n(self):
        if has_equations:
            n = use_equation(self.name, 'n')
        return n

    def set_k(self, absolute_path):
        """ Same as above, but for k
        """
        self.k = make_list_from_file(absolute_path)

    def get_k(self):
        if has_equations:
            k = use_equation(self.name, 'k')
        return k

    def get_has_equations(self):
        return has_equations

    def set_has_equations(self, bool):
        has_equations = bool
    
#     def save_material(self):
#        """Saves the material properties to a file for ease of use in the future"""
#        TO DO

#     def load_material(material_name):
#        """Loads the material properties from a file based on the material_name (string) given"""
#        TO DO

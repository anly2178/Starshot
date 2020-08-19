from .save_load_mat import save_material, del_material, material_exists, load_material, make_list_from_file
from .interpolator import interpolate_from_list
import scipy
import numpy as np
from os import path

""" Each material should have a:
        - name (str; potentially multiple names. Main name is chemical formula
          e.g. SiO2, GeO2, Si3N4, etc.)
        - density (float; kg/m^3)
        - max_temp (float; K)
            * is simply the temperature beyond which a material cannot be
              considered structurally sound. For most materials, this can be
              considered the melting point of the material, but may vary in
              special cases. For example, consider glasses like SiO2, GeO2
              which can become quite viscous past their glass transition
              temperature which lies lower than melting point
        - set of optical constants
            * real component (list of list/tuples, each tuple has 2 parts,
              the wavelength and the value)
            * imaginary component (same as above)
        - (optional) a list of entries detailing equations used for a material
          which detail optical constants
            * each entry is another list of form:
                [ str: {name}
                  list/tuple: {range},
                  func: {equation}
                ]
            * the RANGE in which an equation can be used accurately is described
              by start_wavelength and end_wavelength. Range goes FROM
              start_wavelength TO end_wavelength:
                range = [float {start_wavelength}, float {end_wavelength}]
              BOTH wavelengths should be given in metres.

IMPORTANT NOTE: If using an equation, ensure the equation TAKES IN wavelengths
                IN METRES as units. This is especially important when using
                equations from papers which may use wavenumbers in their equations
                in which case the user is expected to convert the units to
                micrometres within the scope of their defined/loaded function
"""

class Material:

    def __init__(self, name=None, density=None, max_temp=None, abs_coeff=None, n_list_path=None, k_list_path=None):
        """ Constructor requires at least the name and the density
        """
        if material_exists(name):
            mat = load_material(name)
            self.name = mat.get_name()
            self.density = mat.get_density()
            self.max_temp = mat.get_max_temp()
            self.abs_coeff = mat.get_abs_coeff()
            self.n_list = mat.get_n_list()
            self.k_list = mat.get_k_list()
            self.n_equations = mat.get_n_equations()
            self.k_equations = mat.get_k_equations()
        else:
            if name is None:
                raise ValueError('Enter name')
            self.name = name
            if density is None:
                raise ValueError('Enter density')
            self.density = density
            if max_temp is None:
                raise ValueError('Enter max_temp')
            self.max_temp = max_temp
            if abs_coeff is None:
                raise ValueError('Enter abs_coeff')
            self.abs_coeff = abs_coeff
            if n_list_path is None or not path.exists(n_list_path[0]):
                self.n_list = None
            else:
                self.n_list = make_list_from_file(n_list_path)
            if k_list_path is None or not path.exists(k_list_path[0]):
                self.k_list = None
            else:
                self.k_list = make_list_from_file(k_list_path)
            self.n_equations = []
            self.k_equations = []
            save_material(self)

    def get_density(self):
        return self.density

    def set_density(self, density):
        self.density = density
        save_material(self)

    def get_max_temp(self):
        return self.max_temp

    def set_max_temp(self, max_temp):
        self.max_temp = max_temp
        save_material(self)

    def get_abs_coeff(self):
        return self.abs_coeff

    def set_abs_coeff(self, abs_coeff):
        self.abs_coeff = abs_coeff
        save_material(self)

    def get_name(self):
        return self.name

    def set_name(self, name):
        old_name = self.get_name()
        self.name = name
        save_material(self)
        del_material(old_name)

    def set_n_list(self, path_flag):
        """ Takes in a string which is the ABSOLUTE (relative should work too,
            but just to be safe) path of the file that contains the data for n.

            File should be a CSV with each entry organised as:
            [wavelength],[n]

            If any of the entries have a length that is not 2, raises a ValueError.
            WAVELENGTHS ARE IN MICROMETRES (this is the most convenient here,
            and it seems refractiveindex.info tends to record in microns as well
            if any data is taken from there)
        """
        self.n_list = make_list_from_file(path_flag)
        save_material(self)

    def get_n(self, wavelength):
        """ If an equation needs to be used, it will use an equation. Each
            equation in materials_equations is identified by the material name
            (and if there is an equation for n and k, the 'n' and 'k' string
            identifiers in the second argument helps to differentiate).
            Requires wavelength as float to find values
        """
        found_equation = False
        for entry in self.n_equations:
            _, range, equation_func = entry     # unpack entry
            start_wavelength, end_wavelength = range       # unpack range
            # Check if in valid range for equation use
            if wavelength >= start_wavelength and wavelength <= end_wavelength:
                n = equation_func(wavelength)
                found_equation = True
        if not found_equation:
            n = interpolate_from_list(self.n_list, wavelength)
        return n

    def get_n_list(self):
        return self.n_list

    def set_k_list(self, path_flag):
        """ Same as above, but for k
        """
        self.k_list = make_list_from_file(path_flag)
        save_material(self)

    def get_k(self, wavelength):
        found_equation = False
        for entry in self.k_equations:
            _, range, equation_func = entry     # unpack entry
            start_wavelength, end_wavelength = range       # unpack range
            # Check if in valid range for equation use
            if wavelength >= start_wavelength and wavelength <= end_wavelength:
                k = equation_func(wavelength)
                found_equation = True
        if not found_equation:
            k = interpolate_from_list(self.k_list, wavelength)
        return k

    def get_k_list(self):
        return self.k_list

    def get_n_equations(self):
        return self.n_equations

    def get_k_equations(self):
        return self.k_equations

    def add_equation(self, name, range, filepath, n_or_k):
        """ Extracts an equation from a .py file that is valid within a specified
            wavelength range (in m) = [start_wavelength, end_wavelength] and saves it
            as a material attribute. A name should also be attached to the
            equation for the sake of easy removal/editing if required. Appends
            the entry to an equations_list list (n_equations or k_equations
            depending on input).

            Raises an error if invalid name (e.g. repeated names)
            Raises an error if there is an equation with a wavelength range
            which overlaps with the equation being added
        """
        if n_or_k == 'n':
            equations_list = self.n_equations
        elif n_or_k == 'k':
            equations_list = self.k_equations
        # Make a dictionary to store the function in
        dic = {}
        exec(open(filepath).read(), globals(), dic)
        func, = dic.values()    # assign var to function for storage in list entry
        entry = [name, range, func]
        if name not in [eq[0] for eq in equations_list]:
            equations_list.append(entry)
        dic.clear()     # Cleared just in case some memory issues occur
        save_material(self)
        return

    def rmv_equation(self, name, n_or_k):
        """ Removed an equation from a saved equations_list (n_equations or k_equations
            depending on input) based on the name initially given by the user
        """
        if n_or_k == 'n':
            equations_list = self.n_equations
        elif n_or_k == 'k':
            equations_list = self.k_equations
        for entry in equations_list:
            if name == entry[0]:
                equations_list.remove(entry)
        save_material(self)
        return

    def print_variables(self):
        """Print attributes of material
        """
        for var, value in self.__dict__.items():
            print(var, '=', value)
        print('')

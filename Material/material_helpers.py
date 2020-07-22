from Material import Material
import materials_equations

def interpolate_from_list(list, wavelength):
    """ Fills in any values using a linear fit between data points given in
        the files. Also sets the values beyond the intervals given in the list
        to 0.
    """

    # From materials, the given list will be in ascending order, so we can
    # use this to our advantage

    if wavelength > list[-1][0]:
        return 0        # "worst case scenario in terms of temperature"
    elif wavelength < list[0][0]:
        return 0        # also "worst case"
    else:
        # performing linear interpolation between points in the list using the
        # point-gradient formula
        i = 0
        while i < len(list)-1:
            if wavelength == list[i][0]:
                val = list[i][1]
            # the first instance where the wavelength input is larger than a
            # value in the list, we want to stop
            elif wavelength > list[i][0]:
            # when this is true, we know val lies between the i'th and i-1'th entries
                value_interval = (list[i][1], list[i+1][1])
                wavelength_interval = (list[i][0], list[i+1][0])
                m = (val_interval[1]-val_interval[0])/(wavelength_interval[1]-wavelength_interval[0])   # rise over run
                y0 = val_interval[0]
                x0 = wavelength_interval[0]
                val = m*(wavelength - x0) + y0
                return val
            else:
                i += 1

def use_equation(material, n_or_k, wavelength):
    """ For materials with equations, this method will find the equation to use
        based on the material name and whether the equation is used to find n or
        k. Requires the wavelength as well (um). Outputs a value of n or k. Since
        equations may vary greatly depending on sources (might use different
        units) or between materials, each equation should be manually implemented
        by the user if required.
    """
    if n_or_k == 'n':
        if material.name == 'SiO2':
            return materials_equations.SiO2_equations(wavelength).real
        elif material.name == 'GeO2':
            if wavelength <= 5:
                return materials_equations.GeO2_Sellmeier(wavelength)
        elif material.name == 'Si3N4':
            if wavelength <= 1.53846:
                return materials_equations.Si3N4_Sellmeier(wavelength)
        # If the wavelength isn't in the right spot for an equation to be used
        # files must be relied upon
        else:
            return interpolate_from_list(material.n_list, wavelength)

    elif n_or_k == 'k':
        if material.name == 'SiO2':
            return materials_equations.SiO2_equations(wavelength).imag

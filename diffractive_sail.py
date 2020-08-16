#Not very useful yet. This should be improved in the future! 

from Starshot.sail import Sail

class DiffractiveSail(Sail):
    """
    Diffractive lightsails.

    ...

    Attributes
    ----------
    name : str
        A name or code that identifies the sail
    mass : float
        Mass of lightsail (excluding payload) [kg]
    area : float
        Area of lightsail [m^2]
    radius : float
        Radius of lightsail [m]
    s_density : float
        Surface density of lightsail [kg/m^2]
    target : float
        Target speed as fraction of speed of light. E.g. 0.2c
    power : float
        Laser power [W]
    wavelength : float
        Laser wavelength [m]
    W : float
        Figure of merit; the square root of 'Reflectivity-adjusted-area-density'
        as defined by Ilic et al. (2018) [sqrt(g)/m]
    diameter : float
        Diameter of laser array [m]
    angles_coeffs : list of tuples of three floats
        Angle [degrees], reflection efficiency and transmission efficiency of each order.

    Methods (for user)
    ------------------
    def __init__(   self, name=None, mass=None, area=None, target=0.2,
                    wavelength=None, power=100e9, angles_coeffs=None):
        The constructor for DiffractiveSail class
    print_variables()
        Prints the variables of the sail
    change_variables()
        Change the variables of the sail
    calculate_mission()
        Calculates the mission scenario, including distance vs speed vs time.
        A folder is created with 2 txt files and 1 png file.
        1 txt file includes distance, speed and time results, the other txt file
        includes the variables of the mission. The png file includes
        speed vs distance and speed vs time graphs.
    """
    def __init__(   self, name=None, mass=None, area=None, target=0.2,
                    wavelength=None, power=100e9, angles_coeffs=None):
        super().__init__(name=name, mass=mass, area=area, target=target,
            power=power, wavelength=wavelength)
        self.angles_coeffs = angles_coeffs
        self.W = self._find_W()
        self.diameter = self._find_diameter()
        self.print_variables()

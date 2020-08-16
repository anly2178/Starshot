from Starshot.motion import state_vs_t
from Starshot.gaussbeam import find_beam_width, find_frac
from Starshot.results import write_results
import numpy as np
from scipy import integrate

class Sail:
    """
    Lightsails.

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
    reflectance : float
        Absolute reflectance of lightsail
    transmittance : float
        Absolute transmittance of lightsail
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
    __init__(   self, name=None, mass=None, area=None, reflectance=None,
                    target=0.2, power=None, wavelength=1.2e-6)
        The constructor for Sail class
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
    def __init__(   self, name=None, mass=None, area=None, reflectance=None,
                    target=0.2, power=None, wavelength=1.064e-6):
        """The constructor for Sail class

        Parameters
        ----------
        name : str
            A name or code that identifies the sail
        mass : float
            Mass of lightsail (excluding payload) [kg]
        area : float
            Area of lightsail [m^2]
        reflectance : float
            Absolute reflectance of lightsail
        target : float
            Target speed as fraction of speed of light. E.g. 0.2c
        power : float
            Laser power [W]
        wavelength : float
            Laser wavelength [m]

        Returns
        -------
        Sail
            Sail with variables specified by user
        """
        self.name = name #Unique code that identifies the sail
        if name is None:
            raise ValueError("Enter name")
        self.mass = mass #kg
        if mass is None:
            raise ValueError("Enter mass")
        self.area = area #m^2
        if area is None:
            raise ValueError("Enter area")
        self.radius = np.sqrt(area/np.pi) #m
        self.s_density = mass/area #kg/m^2
        self.reflectance = reflectance
        if reflectance is not None:
            self.transmittance = 1 - self.reflectance
        else:
            self.transmittance = None
        self.angles_coeffs = [(0, self.reflectance, self.transmittance)] #degrees
        self.target = target #c
        self.power = power #W
        self.wavelength = wavelength #m
        try:
            self.W = self._find_W() #sqrt(g)/m
            self.diameter = self._find_diameter() #m
        except TypeError:
            self.W = None
            self.diameter = None

    def print_variables(self):
        """Prints the variables of the sail

        Parameters
        ----------
        None required

        Returns
        -------
        None
            Prints variables to output
        """
        for variable, value in self.__dict__.items():
            print(variable, '=', value)
        print('')

    def calculate_mission(self):
        """Calculates the mission scenario, including distance vs speed vs time.
        A folder is created with 2 txt files and 1 png file.
        1 txt file includes distance, speed and time results, the other txt file
        includes the variables of the mission. The png file includes
        speed vs distance and speed vs time graphs.

        Parameters
        ----------
        None required

        Returns
        -------
        None
            Creates folder in current working directory.
        """
        state, time = state_vs_t(self)
        beta, dist = state
        write_results(self, beta, dist, time)

    def _find_fraction(self, dist):
        """Calculates the fraction of laser power incident on the lightsail at
        a distance from the laser array.

        Parameters
        ----------
        Distance from the laser array [m]

        Returns
        -------
        float
            The fraction of laser power incident on the lightsail at the distance.
        """
        radius = self.radius #m
        diameter = self.diameter #m
        wavelength = self.wavelength #m
        beam_width = find_beam_width(diameter, wavelength, dist)
        fraction = find_frac(radius, beam_width, dist)
        return fraction

    def _find_diameter(self):
        """Calculates the diameter of the laser array required to achieve the
        target speed.

        Parameters
        ----------
        None required

        Returns
        -------
        float
            The diameter of the laser array [m].
        """
        c = 2.998e8 #m/s
        wavelength = self.wavelength #m
        mass = self.mass * 1000 #g
        W = self.W #sqrt(g)/m
        power = self.power #W
        diameter = (2*wavelength*c**3*np.sqrt(mass)*W)/(np.sqrt(np.pi)*1000*power) #m
        return diameter

    def _find_effective_R(self):
        """Calculates the effective reflectance.

        Parameters
        ----------
        None required

        Returns
        -------
        float
            Effective reflectance
        """
        eff_R = 0
        for angle, r, t in self.angles_coeffs:
            fac = np.cos(np.deg2rad(angle))
            eff_R += r * fac + t * (1 - fac)
        return eff_R

    def _find_W(self):
        """Calculates the square root of RAAD, W, as defined by Ilic 2018.

        Parameters
        ----------
        None required

        Returns
        -------
        float
            Square root of RAAD, W. [sqrt(g)/m]
        """
        s_density = self.s_density * 1000 #g/m^2
        wavelength = self.wavelength #m
        reflectance = self._find_effective_R()
        #To evaluate the integral we define dW
        def dW(beta, reflectance, s_density, wavelength):
            gamma = 1/np.sqrt(1-beta**2)
            ds_wavelength = wavelength*np.sqrt((1+beta)/(1-beta))
            dW = np.sqrt(s_density)/reflectance * (gamma*beta)/(1-beta)**2
            return dW
        target = self.target
        W, _ = integrate.quad(dW, 0, target, args=(reflectance, s_density, wavelength))
        return W

    # def change_variables(self, **kwargs):
    #     """Changes the variables of the sail, as specified using keyworded args.
    #
    #     Parameters
    #     ----------
    #     Keyworded optional arguments of Attributes as above
    #
    #     Returns
    #     -------
    #     None
    #         Updates variables of sail, as specified by user.
    #     """
    #     allowed_keys = {'name', 'material', 'mass', 'thickness', 'area',
    #                     'radius', 'density', 'abs_coeff', 'absorptance',
    #                     'reflectance', 'transmittance', 'W', 'target', 'max_temp'}
    #     self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)
    #     self.print_variables()

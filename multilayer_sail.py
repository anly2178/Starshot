from Starshot.sail import Sail
from Starshot.optical_constants import n_silica, n_germania
import numpy as np
from tmm.tmm import tmm

class MultilayerSail(Sail):
    """
    Multilayer lightsails.

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
    material : list of str
        List of strings representing the materials in each layer
    thickness : list (of floats)
        Thickness of layers [m]
    max_temp : float
        Maximum temperature of sail [K]
    abs_coeff : float
        Absorption coefficient of lightsail. [cm^-1]
    absorptance : float
        Absolute absorption of lightsail

    Methods (for user)
    ------------------
    def __init__(   name=None, material=None, mass=None, thickness=None,
                    area=None, radius=None, s_density=None, abs_coeff=None,
                    absorptance=None, reflectance=None, transmittance=None,
                    W=None):
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
    def __init__(   self, name=None, material=None, mass=None, thickness=None,
                    area=None, reflectance=None, abs_coeff=None, target=0.2,
                    max_temp=1000, power=None, wavelength=1.2e-6):
        """The constructor for MultilayerSail class

        Parameters
        ----------
        name : str
            A name or code that identifies the sail
        material : list of str
            List of strings representing the materials in each layer
        mass : float
            Mass of lightsail (excluding payload) [kg]
        thickness : list (of floats)
            Thickness of layers [m]
        area : float
            Area of lightsail [m^2]
        reflectance : float
            Absolute reflectance of lightsail
        abs_coeff : float
            Absorption coefficient of lightsail. [cm^-1]
        target : float
            Target speed as fraction of speed of light. E.g. 0.2c
        max_temp : float
            Maximum temperature of sail [K]
        power : float
            Laser power [W]
        wavelength : float
            Laser wavelength [m]

        Returns
        -------
        MultilayerSail
            MultilayerSail with variables specified by user
        """
        super().__init__(name, mass, area, reflectance, target, power, wavelength)
        self.material = material
        if material is None:
            raise ValueError("Enter material(s)")
        self.thickness = thickness #m
        if thickness is None:
            raise ValueError("Enter thickness(es)")
        self.max_temp = max_temp #K
        self.abs_coeff = abs_coeff
        if self.abs_coeff is not None:
            self.absorptance = self._find_absorptance()
            if self.power is None:
                self.power = self._find_power()
        elif self.abs_coeff is None and self.power is None:
            raise ValueError("Enter laser power or absorption coefficient")
        if self.reflectance is None:
            self.reflectance = self._find_reflectance()
        if self.transmittance is None:
            self.transmittance = self._find_transmittance()
        self.W = self._find_W()
        self.diameter = self._find_diameter()
        self.print_variables()

    def _find_structure(self):
        """Creates a list representing the structure of the MultilayerSail.

        Parameters
        ----------
        None required

        Returns
        -------
        list of tuples of two floats
            [(refractive index, -thickness [m]), ...]
        """
        n_s = n_silica(self.wavelength); n_g = n_germania(self.wavelength)
        structure = []
        for material, thickness in zip(self.material, self.thickness):
            if material == 'SiO2':
                structure.append((n_s, -thickness))
            elif material == 'GeO2':
                structure.append((n_g, -thickness))
            elif material == 'gap':
                structure.append((1, -thickness))
        return structure

    def _find_absorptance(self):
        """Calculates absorptance of MultilayerSail.

        Parameters
        ----------
        None required

        Returns
        -------
        float
            Absorptance of MultilayerSail
        """
        #Extract parameters
        wavelength = self.wavelength
        abs_coeff = self.abs_coeff
        #Finds the structure
        structure = self._find_structure()

        im_RI = 1j*wavelength*100*abs_coeff/(4*np.pi)
        temp_struc = [(n+im_RI, t) if n!=1 else (n,t) for (n,t) in structure]

        r_p, t_p, r_s, t_s = tmm(temp_struc, wavelength)
        
        R = ((r_p*np.conj(r_p) + r_s*np.conj(r_s))/2).real
        T = ((t_p*np.conj(t_p) + t_s*np.conj(t_s))/2).real
        A = 1 - R - T
        return A

    def _find_reflectance(self):
        """Calculates reflectance of MultilayerSail.

        Parameters
        ----------
        None required

        Returns
        -------
        float
            Reflectance of MultilayerSail
        """
        #Get parameters
        wavelength = self.wavelength
        target = self.target
        structure = self._find_structure()
        shift = np.sqrt((1+target)/(1-target))
        bandwidth = np.linspace(wavelength, wavelength*shift, 100)
        R_all = []
        for b in bandwidth:
            r_p, _, r_s, _ = tmm(structure, b)
            R_all.append( ((r_p*np.conj(r_p) + r_s*np.conj(r_s))/2) )
        R_avg = (sum(R_all)/100).real
        return R_avg

    def _find_transmittance(params):
        """Calculates transmittance of MultilayerSail.

        Parameters
        ----------
        None required

        Returns
        -------
        float
            Transmittance of MultilayerSail
        """
        #Get parameters
        wavelength = self.wavelength
        target = self.target
        structure = self._find_structure()
        shift = np.sqrt((1+target)/(1-target))
        bandwidth = np.linspace(wavelength, wavelength*shift, 100)
        T_all = []
        for b in bandwidth:
            _, t_p, _, t_s = tmm(structure, b)
            T_all.append( ((t_p*np.conj(t_p) + t_s*np.conj(t_s))/2) )
        T_avg = (sum(T_all)/100).real
        return T_avg

    def _find_power(self):
        """Calculates maximum laser power the MultilayerSail can withstand.

        Parameters
        ----------
        None required

        Returns
        -------
        float
            Maximum power for MultilayerSail [W].
        """
        sb = 5.67e-8 #Stefan-Boltzmann constant

        T_max = self.max_temp #K
        area = self.area
        A = self.absorptance

        T_max = params["max_temp"] #K, max temperature sail can sustain
        A = params["area"] #m^2, area of sail
        absorb = params["absorptance"] #Absolute absorption of sail
        #Assume max temperature occurs at beta=0, dist=0
        beta = 0
        dist = 0
        #Bounds for halving the interval
        P_bb = 2*A*sb*T_max**4 #W, Power incident on black body
        P_high = P_bb/absorb #W, Upper bound; power incident on sail accounting for absorptance
        P_low = P_high/1000 #W, Arbitrary lower bound
        #Create structures to not accidentally destroy original parameters
        params_high = params.copy()
        params_high["power"] = P_high
        params_low = params.copy()
        params_low["power"] = P_low
        #Temperatures at bounds
        T_high = find_one_temp(params_high,beta,dist)
        T_low = find_one_temp(params_low,beta,dist)
        #Midpoint temperature
        T_mid = (T_high+T_low)/2
        #Check the run time
        start_time = time.time()
        while abs(T_mid - T_max) >= 0.0005*T_max: #Maybe the accuracy is too low?
            #If the high temperature is too low for some reason
            if T_high <= T_max:
                params_high["power"] = 2*params_high["power"]
            #If the low temperature is too high for some reason
            if T_low >= T_max:
                params_low["power"] = params_low["power"]/2
            #Midpoint
            P_high = params_high["power"]
            P_low = params_low["power"]
            P_mid = (P_high + P_low)/2
            params_mid = params.copy()
            params_mid["power"] = P_mid
            T_mid = find_one_temp(params_mid,beta,dist)
            if T_mid > T_max:
                params_high["power"] = P_mid
            else:
                params_low["power"] = P_mid
            print(P_mid)
        print("--- %s seconds ---" % (time.time() - start_time))
        P_mid = (params_high["power"] + params_low["power"])/2
        return P_mid

    def _find_temp(self, beta, dist):
        """Calculates temperature of MultilayerSail at a speed and distance.

        Parameters
        ----------
        float
            Beta (v/c)
        float
            Distance away [m]

        Returns
        -------
        float
            Equilibrium temperature of MultilayerSail
        """

        #Extracting the parameters
        material = self.material
        power = self.power #W
        mass = self.mass * 1000 #g
        r = power/mass #W/g
        abs_coeff = self.abs_coeff #cm^-1
        thickness = self.thickness #m
        wavelength_0 = self.wavelength #m
        area = self.area #m^2
        s_density = self.s_density * 1000 #gm^-2
        structure = self._find_structure()
        #Finding power absorbed, accounting for doppler shift and diffraction effects
        wavelength = wavelength_0*np.sqrt((1+beta)/(1-beta))
        A = self.absorptance
        fraction = self._find_fraction(dist=0)
        power_in = r*A*s_density*(1-beta)/(1+beta) * fraction

        """ Note: Honestly, this below section should be made into its own function
            since it is a reusable block of code. Consider doing this at some point
            but for now, focus on optimising code and commenting
        """
        # The RHS is more complicated, since you can't get an expression for T explicitly
        # We need to integrate power flux over all wavelengths to get the total radiated power
        midpoint = 0
        bb_temp = (power_in/(2*1*5.67e-8))**0.25

        T_low = bb_temp             # Lower bound = max emissivity = black body temp
        T_high = bb_temp*10         # Upper bound arbitrary (might not hold at higher temps) - should find a way to set a true reasonable higher bound
        # Use trapezoidal rule to find the total power out for a given temperature
        def power_out(T):
            points = 101            # Can be changed if better resolution is required

            # Ilic paper uses 1-25 microns, but eqns should be valid from 0.5-50 microns if so required
            bounds = np.linspace(1e-6, 25e-6, points)
            power_out_at_wl = np.zeros(points)

            # Running each integral and adding to the list (optimisation here would be to fix list size and assign vals)
            i = 0
            for wavelength in bounds:
                poawl = (spectral_power_flux(wavelength, structure, T))
                if poawl == None: #In the case of overflow error, return None
                    return None
                else:
                    power_out_at_wl[i] = poawl
                i += 1
            power_out = np.trapz(power_out_at_wl, bounds, (25e-6-1e-6)/points)
            return power_out

        start_time = time.time()
        # Powers at the bounds of temperature interval
        P_high = power_out(T_high)
        P_low = power_out(T_low)
        if P_high == None or P_low == None: #If overflow error, return temperature as 0
            return 0

        # Halving the interval for a result
        while abs(P_high - P_low) >= 0.01*power_in:
        # The only issue we can really get is if P_high is too low - if this is
        # the case, just double P_high
            if (P_high <= power_in):
                T_high = T_high*2

            midpoint = (T_low+T_high)/2
            P_mid = power_out(midpoint)
            if P_mid == None: #If overflow error, return temperature as 0
                return 0
            if P_mid > power_in:
                T_high = midpoint
            else:
                T_low = midpoint

            P_high = power_out(T_high)
            P_low = power_out(T_low)
            if P_high == None or P_low == None: #If overflow error, return temperature as 0
                return 0
            print(midpoint)
        # Take the midpoints as the final result since this is the result from halving interval
        midpoint = (T_high+T_low)/2
        print("--- %s seconds ---" % (time.time() - start_time))
        return midpoint
    #Methods here

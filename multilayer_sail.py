from Starshot.sail import Sail
from Starshot.tmm.tmm import tmm
from Starshot.materials.save_load_mat import load_material
import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi

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
    materials : list of str
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
    def __init__(   name=None, materials=None, mass=None, thickness=None,
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
    def __init__(   self, name=None, materials=None, mass=None, thickness=None,
                    area=None, reflectance=None, abs_coeff=None, target=0.2,
                    max_Starchip_temp=1000, power=None, wavelength=1.064e-6):
        """The constructor for MultilayerSail class
        Parameters
        ----------
        name : str
            A name or code that identifies the sail
        materials : list of str
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
        max_Starchip_temp : float
            Maximum temperature of sail payload [K]
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
        # This block converts the material names (list of strings) into a list of
        # Material objects
        if materials is None:
            raise ValueError("Enter material(s)")
        try:
            self.materials = [load_material(mat) for mat in materials]
        except ValueError as ve:
            print(ve)
        self.thickness = thickness #m
        if thickness is None:
            raise ValueError("Enter thickness(es)")
        self.max_Starchip_temp = max_Starchip_temp #K
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

    def _find_structure(self, wavelength = self.wavelength):
        """Creates a list representing the structure of the MultilayerSail.
        Parameters
        ----------
        None required
        Returns
        -------
        list of tuples of two floats
            [(refractive index, -thickness [m]), ...]
        """
        structure = []
        for material, thickness in zip(self.material, self.thickness):
            structure.append( (material.get_n(wavelength) + 1j*material.get_k(wavelength), -thickness) )
        return structure

    def _find_SA_density(self):
        """ Determines the surface area density of the sail given its structure
            and mass densities of each material
            Parameters
            ----------
            None required
            Returns
            ----------
            float
                surface area density [kg m-2]
        """
        structure = self.structure
        SA_density = 0
        for material, thickness in zip(self.materials, self.thickness):
            SA_density += material.get_density()*thickness
        return SA_density

    def _find_absorptance(self, wavelength = self.wavelength):
        """Calculates absorptance of MultilayerSail based on the (expected)
        absorption coefficients of the sail materials (material.abs_coeff
        attribute) and a wavelength being analysed within the laser bandwidth.
        Assumes the absorption coefficient is constant over the laser range
        in the near IR. Usage of absorption coefficient in the laser bandwidth
        is due to extinction coefficients for most materials not being well
        established in this range. Further, results for extinction and
        absorption coefficients may vary depending on purity and manufacture of
        materials.
        Parameters
        ----------
        float (optional)
            wavelength [m]
        Returns
        -------
        float
            Absorptance of MultilayerSail
        """
        structure_near_IR = []
        for material, thickness in zip(self.materials, self.thickness):
            k = 1j*wavelength*100*material.get_abs_coeff()/(4*pi)   # conversion from abs_coeff to extinction coeff
            structure_near_IR.append( (material.get_n(wavelength) + k, -thickness) )

        r_p, t_p, r_s, t_s = tmm(structure_near_IR, wavelength)

        R = ((r_p*np.conj(r_p) + r_s*np.conj(r_s))/2).real
        T = ((t_p*np.conj(t_p) + t_s*np.conj(t_s))/2).real
        A = 1 - R - T
        return A

    def _find_reflectance(self):
        """Calculates reflectance of MultilayerSail, averaged over wavelength.
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
        """Calculates transmittance of MultilayerSail, averaged over wavelength.
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

    def _spectral_power_flux(self, wavelength, temperature, points_in_integration = 50):

        """ Finds the spectral power flux of an "ideal" (perfectly flat and smooth)
        sail that is SYMMETRICAL on both sides of the sail. This is the energy
        emitted per unit area at given wavelength. Accounts for asymmetric
        multilayer_sails (along the axis of the incident laser light)

        Uses trapezoidal rule to integrate power emitted at all angles to
        produce the spectral power flux

        Parameters
        ----------
        float
            wavelength [m]
        float
            temperature [m]
        int
            points_in_integration
                - this is the number of points used in the trapezoidal rule
                  integration
        Returns
        -------
        float
            emissivity in direction described by angle and sail structure
        """

        def directional_emissivity(self, angle, wavelength, front_or_back):
            """ Calculates the directional emissivity of a given multilayer_sail
                structure based on a wavelength (float) and incident angle (i.e.
                angle of elevation). This assumes the sail is perfectly smooth and
                the structure is radially symmetric along the surface of the sail at
                each point of the sail.
                Parameters
                ----------
                float
                    angle [radians]
                float
                    wavelength [m]
                string
                    front_or_back
                        - denotes whether we are calculating the directional
                          emissivity of the front or back face of the sail
                Returns
                -------
                float
                    emissivity in direction described by angle and sail structure
            """
            # Creates a new structure list that is based on calculated optical constants using wavelength
            if front_or_back == 'front':
                structure = self._find_structure(wavelength)
            elif front_or_back == 'back':
                structure = self._find_structure(wavelength).reverse()

            # This block gives an expression for emissivity in terms of theta (and wavelength)
            # First set out by finding the reflectance and transmittance of structure at
            # this wavelength and angle
            r_p, t_p, r_s, t_s = tmm(structure, wavelength, theta)
            R = ( r_p*np.conj(r_p) + r_s*np.conj(r_s) )/2
            T = ( t_p*np.conj(t_p) + t_s*np.conj(t_s) )/2
            dEpsilon = (1-R-T)
            return dEpsilon.real

        # First, give expression for radiation emitted by a black body
        h = 6.62607004e-34       # Planck's constant in SI
        c = 299792458             # speed of light in SI
        k_B = 1.38064852e-23        # Boltzmann constant in SI
        I = ((2*h*c**2)/wavelength**5)*(1/(np.exp(h*c/(wavelength*k_B*temperature))-1))         # Planck's Law

        # Now give expression for hemispherical emissivity. Note factor of 2: 2 comes
        # from integrating wrt phi (the azimuth)

        # Use trapezoidal integration to speed things up
        bounds = np.linspace(0,pi/2,points_in_integration)

        # ONCE FOR EMISSION FROM FRONT FACE
        direc_ems = [2*directional_emissivity(theta, wavelength, 'front')*cos(theta)*sin(theta) for theta in bounds]
        # In the below line, note that the integration returns the spectral hemispherical emissivity
        front_power_flux = pi*I*np.trapz(direc_ems, bounds, pi/2/points_in_integration)

        # SECOND TIME FOR BACK FACE
        direc_ems = [2*directional_emissivity(theta, wavelength, 'back')*cos(theta)*sin(theta) for theta in bounds]
        back_power_flux = pi*I*np.trapz(direc_ems, bounds, pi/2/points_in_integration)

        power_flux = front_power_flux + back_power_flux
        return power_flux

    def _find_eq_temps_given_abs_coeff(self):
        """ Determines the maximum equilibrium temperature of the sail given
            the absorption coefficients of each material in the sail.
            If any of the materials do not have an allocated absorption
            coefficient, will raise an exception.
            Parameters
            ----------
            None required
            Returns
            -------
            float
                equilibrium temperature [K]
        """
        initial_wavelength = self.wavelength        # laser wavelength
        target = self.target
        structure = self.structure.copy()
        # Below block of code finds the maximum power absorbed by the sail throughout its journey
        betas = np.linspace(0,target,100)  # fraction of speed of light sail is travelling at
        power_absorbed = 0       # need to find maximum p_in based on beta
        power_mass_ratio = self.power/self.mass     # laser power to mass ratio of sail
        s_density = self.s_density         # surface area density

        # Loop that gets the maximum power value
        for beta in betas:
            wavelength = initial_wavelength*np.sqrt((1+beta)/(1-beta))
            A = self._find_absorptance(wavelength)

            # Finding the LHS of Atwater et al. 2018's  equation
            power_beta = power_mass_ratio*A*s_density*(1-beta)/(1+beta)       # power absorbed when v/c = beta

            if power_beta > power_absorbed:
                power_absorbed = power_beta     # since maximum power in results in highest equilibrium temperature

        def power_in_minus_out(T, power_absorbed):
            """ Uses an input temperature to find the total power emitted by the
                sail. Subtracts this value from the power absorbed, given as input.
                Roots occur when the power in = power out, and hence at thermal
                equilibrium.
                Parameters
                ----------
                float
                    T(emperature) [K]
                float
                    power_absorbed []
                Returns
                ----------
                float
                    difference []
                        - difference between power_absorbed and power_emitted
            """

            def find_power_emitted(T, points_in_integration = 100, integration_range = [1e-6, 25e-6]):
                """ Finds the power emitted by a sail with given structure at a
                    specific temperature. Determnied by performing a trapezoidal
                    integration over a (default 1-25 micron) wavelength range of the
                    spectral power flux, calculated by the _spectral_power_flux()
                    method.
                    Parameters
                    ----------
                    float
                        T(emperature) [K]
                    int (optional)
                        points_in_integration
                            - number of points used in trapezoidal integration
                    list/tuple (optional)
                        integration_range
                            - wavelength range over which spectral power flux
                              is integrated over to determine total power per
                              unit area of sail emitted (note the area of the
                              sail in this respect is the area of one face,
                              NOT the surface area = 2 * sail area)
                    Returns
                    ----------
                    float
                        power_emitted []
                """
                lower_bound, upper_bound = integration_range
                points = np.linspace(lower_bound, upper_bound, points_in_integration)
                # Calling _spectral_power_flux at each point and adding to the list for integration
                power_out_at_wl = [self._spectral_power_flux(wavelength,T) for wavelength in points]
                power_emitted = np.trapz(power_out_at_wl, points, (upper_bound - lower_bound)/points_in_integration)
                return power_emitted

            return power_absorbed - find_power_emitted(T)

        # The zero of the _power_in_minus_out function occurs when the sail is at
        # equilibrium temperature. Hence, we can use Newton's method to estimate
        # this temperature

        # Use the starting point of Newton's method as the black body temperature
        # given the power_absorbed (per unit area)

        bb_temp = (power_in/(2*1*5.67e-8))**0.25

        eq_temp = scipy.optimize.newton(power_in_minus_out, bb_temp, args = (power_absorbed))

        return eq_temp

# ============================================================================
# Going to cut off these functions here for now - I'm not too sure what they do
# and they're a bit technical, so I'll leave these untouched for now, sorry.

#I will comment these out for now so they do not interfere. If we find them
#unnecessary in the future, we will delete them.
    #
    # def _find_temp(self, beta, dist):
    #     """Calculates temperature of MultilayerSail at a speed and distance.
    #     Parameters
    #     ----------
    #     float
    #         Beta (v/c)
    #     float
    #         Distance away [m]
    #     Returns
    #     -------
    #     float
    #         Equilibrium temperature of MultilayerSail
    #     """
    #
    #     #Extracting the parameters
    #     material = self.material
    #     power = self.power #W
    #     mass = self.mass * 1000 #g
    #     r = power/mass #W/g
    #     abs_coeff = self.abs_coeff #cm^-1
    #     thickness = self.thickness #m
    #     wavelength_0 = self.wavelength #m
    #     area = self.area #m^2
    #     s_density = self.s_density * 1000 #gm^-2
    #     structure = self._find_structure()
    #     #Finding power absorbed, accounting for doppler shift and diffraction effects
    #     wavelength = wavelength_0*np.sqrt((1+beta)/(1-beta))
    #     A = self.absorptance
    #     fraction = self._find_fraction(dist=0)
    #     power_in = r*A*s_density*(1-beta)/(1+beta) * fraction
    #
    #     """ Note: Honestly, this below section should be made into its own function
    #         since it is a reusable block of code. Consider doing this at some point
    #         but for now, focus on optimising code and commenting
    #     """
    #     # The RHS is more complicated, since you can't get an expression for T explicitly
    #     # We need to integrate power flux over all wavelengths to get the total radiated power
    #     midpoint = 0
    #     bb_temp = (power_in/(2*1*5.67e-8))**0.25
    #
    #     T_low = bb_temp             # Lower bound = max emissivity = black body temp
    #     T_high = bb_temp*10         # Upper bound arbitrary (might not hold at higher temps) - should find a way to set a true reasonable higher bound
    #     # Use trapezoidal rule to find the total power out for a given temperature
    #     def power_out(T):
    #         points = 101            # Can be changed if better resolution is required
    #
    #         # Ilic paper uses 1-25 microns, but eqns should be valid from 0.5-50 microns if so required
    #         bounds = np.linspace(1e-6, 25e-6, points)
    #         power_out_at_wl = np.zeros(points)
    #
    #         # Running each integral and adding to the list (optimisation here would be to fix list size and assign vals)
    #         i = 0
    #         for wavelength in bounds:
    #             poawl = (spectral_power_flux(wavelength, structure, T))
    #             if poawl == None: #In the case of overflow error, return None
    #                 return None
    #             else:
    #                 power_out_at_wl[i] = poawl
    #             i += 1
    #         power_out = np.trapz(power_out_at_wl, bounds, (25e-6-1e-6)/points)
    #         return power_out
    #
    #     start_time = time.time()
    #     # Powers at the bounds of temperature interval
    #     P_high = power_out(T_high)
    #     P_low = power_out(T_low)
    #     if P_high == None or P_low == None: #If overflow error, return temperature as 0
    #         return 0
    #
    #     # Halving the interval for a result
    #     while abs(P_high - P_low) >= 0.01*power_in:
    #     # The only issue we can really get is if P_high is too low - if this is
    #     # the case, just double P_high
    #         if (P_high <= power_in):
    #             T_high = T_high*2
    #
    #         midpoint = (T_low+T_high)/2
    #         P_mid = power_out(midpoint)
    #         if P_mid == None: #If overflow error, return temperature as 0
    #             return 0
    #         if P_mid > power_in:
    #             T_high = midpoint
    #         else:
    #             T_low = midpoint
    #
    #         P_high = power_out(T_high)
    #         P_low = power_out(T_low)
    #         if P_high == None or P_low == None: #If overflow error, return temperature as 0
    #             return 0
    #         print(midpoint)
    #     # Take the midpoints as the final result since this is the result from halving interval
    #     midpoint = (T_high+T_low)/2
    #     print("--- %s seconds ---" % (time.time() - start_time))
    #     return midpoint
    #Methods here


# Not sure what this is meant to do, moving to bottom of file

    # def _find_power(self):
    #     """Calculates maximum laser power the MultilayerSail can withstand.
    #     Parameters
    #     ----------
    #     None required
    #     Returns
    #     -------
    #     float
    #         Maximum power for MultilayerSail [W].
    #     """
    #     sb = 5.67e-8 #Stefan-Boltzmann constant
    #
    #     T_max = self.max_temp #K
    #     area = self.area
    #     A = self.absorptance
    #
    #     T_max = params["max_temp"] #K, max temperature sail can sustain
    #     A = params["area"] #m^2, area of sail
    #     absorb = params["absorptance"] #Absolute absorption of sail
    #     #Assume max temperature occurs at beta=0, dist=0
    #     beta = 0
    #     dist = 0
    #     #Bounds for halving the interval
    #     P_bb = 2*A*sb*T_max**4 #W, Power incident on black body
    #     P_high = P_bb/absorb #W, Upper bound; power incident on sail accounting for absorptance
    #     P_low = P_high/1000 #W, Arbitrary lower bound
    #     #Create structures to not accidentally destroy original parameters
    #     params_high = params.copy()
    #     params_high["power"] = P_high
    #     params_low = params.copy()
    #     params_low["power"] = P_low
    #     #Temperatures at bounds
    #     T_high = find_one_temp(params_high,beta,dist)
    #     T_low = find_one_temp(params_low,beta,dist)
    #     #Midpoint temperature
    #     T_mid = (T_high+T_low)/2
    #     #Check the run time
    #     start_time = time.time()
    #     while abs(T_mid - T_max) >= 0.0005*T_max: #Maybe the accuracy is too low?
    #         #If the high temperature is too low for some reason
    #         if T_high <= T_max:
    #             params_high["power"] = 2*params_high["power"]
    #         #If the low temperature is too high for some reason
    #         if T_low >= T_max:
    #             params_low["power"] = params_low["power"]/2
    #         #Midpoint
    #         P_high = params_high["power"]
    #         P_low = params_low["power"]
    #         P_mid = (P_high + P_low)/2
    #         params_mid = params.copy()
    #         params_mid["power"] = P_mid
    #         T_mid = find_one_temp(params_mid,beta,dist)
    #         if T_mid > T_max:
    #             params_high["power"] = P_mid
    #         else:
    #             params_low["power"] = P_mid
    #         print(P_mid)
    #     print("--- %s seconds ---" % (time.time() - start_time))
    #     P_mid = (params_high["power"] + params_low["power"])/2
    #     return P_mid

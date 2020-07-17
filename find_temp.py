def find_one_temp(sail, beta, dist):
    """
    Returns the equilibrium temperature at a specific speed and distance,
    accounting for doppler shift and diffraction effects.

    Input:  parameters, beta (fraction of c), distance from DE system (m)
    Ouput:  Equilibrium temperature (K)
    """

    #Extracting the parameters
    material = sail.material
    power = sail.power #W
    mass = sail.mass * 1000 #g
    r = power/mass #W/g
    abs_coeff = sail.abs_coeff #cm^-1
    thickness = sail.thickness #m
    wavelength_0 = sail.wavelength #m
    area = sail.area #m^2
    s_density = sail.s_density * 1000 #gm^-2

    #Determining which material to calculate absorption for
    if material == 'SiO2':
        n = 1.45
    elif material == 'GeO2':
        n = 1.6
    #Create structure
    i = 0
    structure = []
    while i < len(thickness):
        if i%2 == 0:
            layer = (n,-thickness[i]) #Material
        elif i%2 == 1:
            layer = (1,-thickness[i]) #Vacuum
        structure.append(layer)
        i += 1
    #Finding power absorbed, accounting for doppler shift and diffraction effects
    wavelength = wavelength_0*np.sqrt((1+beta)/(1-beta))
    A = find_absorption_from_coefficient(structure, abs_coeff, wavelength)
    try:
        fraction = find_fraction_incident(params, dist)
    except TypeError: #Diameter not given, which occurs when finding max power
        fraction = 1 #Very good approx for small distances
    power_in = ratio*A*rho_S*(1-beta)/(1+beta) * fraction

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

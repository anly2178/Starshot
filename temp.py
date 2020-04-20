from .TMM_analysis_sail.find_eq_temp import *
from .TMM_analysis_sail.optical_constants import n_silica

def solveTemp(params, abs_coeff, beta):
    # Absorbance is measured as 1-R-T where R, T are reflectance and transmittance.
    # dependence on the absorption coefficient is based on definition of RI from paper.

    """ TURNS OUT, WE HAVE TO ACCOUNT FOR THE DOPPLER SHIFT!
    THE ILIC PAPER PLOTS THE MAXIMUM TEMPERATURE CURVE BASED ON THE SPEED OF
    THE CRAFT. HENCE WE HAVE TO FIND THE VALUE FOR BETA WHERE THE POWER IN PER
    UNIT AREA IS MAXIMISED AS WELL. SO WE RUN THIS THROUGH CODE BELOW:
    """

    laser_power = params["power"] #W
    m_sail = params["m_sail"] * 1e3 #g
    ratio = laser_power/m_sail #W/g

    thickness = params["thickness"] #m
    wavelength_0 = params["wavelength"] #m
    n = n_silica(wavelength_0)
    structure = [(n,-thickness)]

    density = params["density"]
    rho_S = density * thickness

    #Finding power absorbed, accouting for doppler shift
    wavelength = wavelength_0*np.sqrt((1+beta)/(1-beta))
    A = find_absorption_from_coefficient(structure, abs_coeff, wavelength)
    power_in = ratio*A*rho_S*(1-beta)/(1+beta)

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
            power_out_at_wl[i] = (spectral_power_flux(wavelength, structure, T))
            i += 1
        power_out = np.trapz(power_out_at_wl, bounds, (25e-6-1e-6)/points)
        return power_out

    start_time = time.time()
    # Powers at the bounds of temperature interval
    P_high = power_out(T_high)
    P_low = power_out(T_low)

    # Halving the interval for a result
    while abs(P_high - P_low) >= 0.05*power_in:

    # The only issue we can really get is if P_high is too low - if this is
    # the case, just double P_high
        if (P_high <= power_in):
            T_high = T_high*2

        midpoint = (T_low+T_high)/2

        if power_out(midpoint) > power_in:
            T_high = midpoint
        else:
            T_low = midpoint

        P_high = power_out(T_high)
        P_low = power_out(T_low)
    # Take the midpoints as the final result since this is the result from halving interval
    midpoint = (T_high+T_low)/2
    print("--- %s seconds ---" % (time.time() - start_time))
    return midpoint

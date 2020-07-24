import scipy
import numpy as np

""" Any equations that are required to describe the optical constants of a given
    material should be included here (when a table of values cannot be used).
    Existing implementations of equations are used here as examples. If these
    equations are only useful in a specific wavelength range, see the
    material_helpers.use_equation() method for further details on implementation
"""

def SiO2_equations(wavelength):

    """ Returns refractive index of SiO_2 from wavelength argument.
    Input wavelength is in metres. Sellmeier equation is only valid
    from 1-8 micrometres, but that is good enough for calculating the
    reflectivity of a silica-based lightsail for the laser bandwidth.
    Source: Rei Kitamura, Laurent Pilon and Miroslaw Jonasz (2007).

    Past Sellmeier range, use Kitamura, Pilon, Jonasz's eqn for 7-50 um
    """

    """ Nested function defining some functions required to find the complex
        refractive index of silica. Wavenumber is in cm^-1 """

    def silica_g(wavenumber, is_kkg = False):

        # Quickly defining another function that is requred later on
        def D(x):
            return np.array([scipy.special.dawsn(y) for y in x])

        alpha = np.array((3.7998, 0.46089, 1.2520, 7.8147, 1.0313, 5.3757, 6.3305, 1.2948))
        wavenumber0 = np.array((1089.7, 1187.7, 797.78, 1058.2, 446.13, 443, 465.80, 1026.7))
        sigma = np.array((31.454, 100.46, 91.601, 63.153, 275.111, 45.220, 22.680, 232.14))

        if is_kkg == True:
            g = (2/pi)*alpha*( D( 2 * np.sqrt(np.log(2)) * ((wavenumber+wavenumber0)/sigma) ) - D( 2 * np.sqrt(np.log(2)) * ((wavenumber-wavenumber0)/sigma) ) )
        else:
            g = alpha*(np.exp(-4*np.log(2)*((wavenumber-wavenumber0)/sigma)**2) - np.exp(-4*np.log(2)*((wavenumber+wavenumber0)/sigma)**2))
        return g

    # Assume k negligible below 7 microns; use Sellmeier equation

    if wavelength < 7e-6:
        sellmeier_wl = (wavelength*1e6)**2
        term_a = (0.6961663*sellmeier_wl)/(sellmeier_wl-(0.0684043)**2)
        term_b = (0.4079426*sellmeier_wl)/(sellmeier_wl-0.1162414**2)
        term_c = (0.8974794*sellmeier_wl)/(sellmeier_wl-9.896161**2)
        n = np.sqrt(1 + term_a + term_b + term_c)
    # k no longer negligible
    elif wavelength >= 7e-6:
        wavenumber = (1/wavelength)/100     # in cm^-1
        g_kkg = silica_g(wavenumber, True)
        g = silica_g(wavenumber)
        n = np.sqrt(2.1232 + np.sum(g_kkg) + 1j*np.sum(g))
    return n

def GeO2_Sellmeier(wavelength):

    """ Provides a Sellmeier equation for GeO2 which is accurate below 5 microns
    Source: Fleming, J. W., 1984
    """

    A1 = 0.80686642
    A2 = 0.71815848
    A3 = 0.85416831
    l1 = 0.68972606e-1
    l2 = 0.15396605
    l3 = 0.11841931e2
    wlsq = wavelength**2
    nsq = 1 + (wlsq*A1)/(wlsq-l1**2) + (wlsq*A2)/(wlsq-l2**2) + (wlsq*A3)/(wlsq-l3**2)
    return np.sqrt(nsq)

def Si3N4_Sellmeier(wavelength):

    """ Provides a Sellmeier equation for Si3N4 which is accurate below ~1.5 microns
    Source: Luke et al., 2015
    """

    wavelength = wavelength
    A1 = 3.0249
    A2 = 40314
    l1 = 0.1353406
    l2 = 1239.842
    wlsq = (wavelength)**2
    nsq = 1 + (wlsq*A1)/(wlsq-l1**2) + (wlsq*A2)/(wlsq-l2**2)
    return np.sqrt(nsq)

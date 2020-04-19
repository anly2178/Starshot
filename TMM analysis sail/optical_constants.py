import scipy
import scipy.integrate as integrate
import numpy as np
from numpy import sin, cos, pi

""" Returns refractive index of SiO_2 from wavelength argument.
    Input wavelength is in metres. Sellmeier equation is only valid
    from 1-8 micrometres, but that is good enough for calculating the
    reflectivity of a silica-based lightsail for the laser bandwidth.
    Source: Rei Kitamura, Laurent Pilon and Miroslaw Jonasz (2007).

    Past Sellmeier range, use Kitamura, Pilon, Jonasz's eqn for 7-50 um
"""
def n_silica(wavelength):

    """ Nested function defining some functions required to find the complex
        refractive index of silica. Wavenumber is in cm^-1 """

    def silica_g(wavenumber, is_kkg = False):

        # Quickly defining another function that is requred later on
        def D(x):
            D = []
            i = 0
            while i < len(x):
                d = scipy.special.dawsn(x[i])
                D.append(d)
                i += 1
            return np.array(D)

        alpha = np.array((3.7998, 0.46089, 1.2520, 7.8147, 1.0313, 5.3757, 6.3305, 1.2948))
        wavenumber0 = np.array((1089.7, 1187.7, 797.78, 1058.2, 446.13, 443, 465.80, 1026.7))
        sigma = np.array((31.454, 100.46, 91.601, 63.153, 275.111, 45.220, 22.680, 232.14))
        if is_kkg == True:
            g = (2/pi)*alpha*( D( 2 * np.sqrt(np.log(2)) * ((wavenumber+wavenumber0)/sigma) ) - D( 2 * np.sqrt(np.log(2)) * ((wavenumber-wavenumber0)/sigma) ) )
        else:
            g = alpha*(np.exp(-4*np.log(2)*((wavenumber-wavenumber0)/sigma)**2) - np.exp(-4*np.log(2)*((wavenumber+wavenumber0)/sigma)**2))
        return g

    """ Back to the main function """
    # Assume k negligible
    if wavelength < 7e-6:
        sellmeier_wl = (wavelength*1e6)**2
        term_a = (0.6961663*sellmeier_wl)/(sellmeier_wl-(0.0684043)**2)
        term_b = (0.4079426*sellmeier_wl)/(sellmeier_wl-0.1162414**2)
        term_c = (0.8974794*sellmeier_wl)/(sellmeier_wl-9.896161**2)
        n = np.sqrt(1 + term_a + term_b + term_c)
    # k no longer negligible
    elif wavelength >= 7e-6:
        wavenumber = (1/wavelength)/100     # in cm^-1
        sum = 0
        g_kkg = silica_g(wavenumber, True)
        g = silica_g(wavenumber)
        n = np.sqrt(2.1232 + np.sum(g_kkg) + 1j*np.sum(g))
    return n


""" Calculating the refractive index of crystalline silicon given different
    wavelengths using a 3 term Sellmeier equation.
    Source: Horowitz, Amirtharaj (2004)

    *ACCEPTS NO ARGUMENT* This is because only discrete datapoints can be
    obtained and the result cannot be determined by an equation (that I have found)

    Output: List of size 294 (since this divides data points well), which can be
    used for analysis. Difficult to work with due to divisibility of the number
    of points given, but code can be written to work around this. List contains
    tuples of kind (wavelength (in MICRONS), complex refractive index)
"""

def n_csi():
    # Equation and constants given are for wavenumbers (cm^-1), so we need to convert
    def get_n(wl):
        wl_in_cm = wl*1e-4
        wavenumber = 1/wl_in_cm
        wnsq = wavenumber**2
        # Sellmeier constants
        eps = 11.67316
        A = 10e-9
        B = 0.00365
        C = 9.0236e3
        nsq = eps + A*wnsq + (B*wnsq)/(C**2-wnsq)
        return np.sqrt(nsq)

    # For k, read text file, take n elements and return a list of n elements
    # linearly spaced
    def get_list_k(n):
        f = open('k_silicon.txt', 'r')
        points = f.readlines()
        interval = len(points)/(n-1)
        i = 0
        k = []
        while i < len(points):
            string = points[int(i)].strip()
            point = string.split(",,")
            k.append((float(point[0]),float(point[1])))
            i+=interval
        return k

    k_list = (get_list_k(294))

    n = []
    for point in k_list:
        n.append((point[0],get_n(point[0])+point[1]*1j))

    return n

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

""" Calculating the complex refractive index of GeO2 at different
    wavelengths. Valid from 1-50 um (with interpolation between 4.8-5 um
    using Fleming's Sellmeier equation.
"""

def n_germania(wavelength):
    # Working in um
    wavelength = wavelength*1e6
    # Functions for reading files
    def get_list_n():
        f = open('starshot_lib/TMM_analysis_sail/n_germania.txt', 'r')
        points = f.readlines()
        i = 0
        n = []
        while i < len(points):
            string = points[i].strip()
            point = string.split('\t')
            n.append((float(point[0]), float(point[1])))
            i += 1
        f.close()
        return n

    def get_list_k():
        f = open('starshot_lib/TMM_analysis_sail/k_germania.txt', 'r')
        points = f.readlines()
        i = 0
        k = []
        while i < len(points):
            string = points[i].strip()
            point = string.split('\t')
            k.append((float(point[0]), float(point[1])))
            i += 1
        f.close()
        return k

    # Now to actually get what we need to return
    def find_val_from_list(ls, wavelength):
        # Each list is in descending order in terms of wavelength, so use this
        # to our advantage. First, we need to find where the wavelength we want
        # lies in the list:
        i = 0
        while i < len(ls):
            if wavelength < ls[-1][0]:
                return ls[-1][1]
            if wavelength > ls[i][0]:
                # Shouldn't be looking past the wavelength ranges that we have here anyways
                if wavelength > ls[0][0]:
                    return 0
                else:
                    val_interval = (ls[i][1],ls[i-1][1])   # (val at lower wl, val at higher wl)
                    wl_interval = (ls[i][0],ls[i-1][0])
                    m = (val_interval[1]-val_interval[0])/(wl_interval[1]-wl_interval[0])
                    y0 = val_interval[0]
                    x0 = wl_interval[0]
                    val = m*(wavelength - x0) + y0
                    return val
            else:
                i += 1

    # And for the range where the above approximation does not hold, and we must
    # instead use a Sellmeier relationship
    def sellmeier(wavelength):
        A1 = 0.80686642
        A2 = 0.71815848
        A3 = 0.85416831
        l1 = 0.68972606e-1
        l2 = 0.15396605
        l3 = 0.11841931e2
        wlsq = (wavelength)**2
        nsq = 1 + (wlsq*A1)/(wlsq-l1**2) + (wlsq*A2)/(wlsq-l2**2) + (wlsq*A3)/(wlsq-l3**2)
        return np.sqrt(nsq)

    if wavelength > 5:
        ns = get_list_n()
        ks = get_list_k()
        return find_val_from_list(ns, wavelength) + 1j*find_val_from_list(ks, wavelength)
    # Need to implement for 1-5 um regime
    elif wavelength <= 5:
        return sellmeier(wavelength)
    else:
        return 0

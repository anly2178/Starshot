def n_silica(wavelength):
    """ This function describes the real refractive index of silica from 1 to 50 micrometres.
        Source: Rei Kitamura, Laurent Pilon and Miroslaw Jonasz (2007).
        doi: 10.1364/AO.46.008118
    NOTE: This function actually describes two separate equations, but includes
              functionality to switch between either function depending on the
	      wavelength given as argument. Such implementations are acceptable for
	      use in this library
    """

    # NOTE: All imports are written WITHIN the scope of the function

    import scipy
    import scipy.integrate as integrate
    import numpy as np
    from numpy import sin, cos, pi

    def silica_g(wavenumber, is_kkg = False):
        """Nested function defining some functions required to find the complex refractive index of silica.
        Wavenumber is in cm^-1 """
        def D(x):
            D = [None]*len(x)
            i = 0
            while i < len(x):
                d = scipy.special.dawsn(x[i])
                D[i] = (d)
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

    if wavelength < 7e-6:
        sellmeier_wl = (wavelength*1e6)**2	# NOTE: Conversion from metres to micrometres (as required for Sellmeier equation abstracted from literature)
        term_a = (0.6961663*sellmeier_wl)/(sellmeier_wl-(0.0684043)**2)
        term_b = (0.4079426*sellmeier_wl)/(sellmeier_wl-0.1162414**2)
        term_c = (0.8974794*sellmeier_wl)/(sellmeier_wl-9.896161**2)
        n = np.sqrt(1 + term_a + term_b + term_c)

    elif wavelength >= 7e-6:
        wavenumber = (1/wavelength)/100     # NOTE: Necessary conversion from wavelength (m) to wavenumber (cm^-1)
        sum = 0
        g_kkg = silica_g(wavenumber, True)
        g = silica_g(wavenumber)
        n = np.sqrt(2.1232 + np.sum(g_kkg) + 1j*np.sum(g))

    return n.real	# taking real component as result since this file should describe the real refractive index only

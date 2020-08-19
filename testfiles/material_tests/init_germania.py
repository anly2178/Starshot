from Starshot.materials.material import Material

# #Initialise germania. Saves GeO2.pkl into saved_materials directory. Note that path is given to a file
# #containing refractive index/extinction coefficient vs wavelength. A flag is also given (1, 2 or 3).
# #Flag 2 is given since wavelengths are in microns.
# #1 is for metres, 2 is for nanometres, 4 is for wavenumber.
germania = Material(name="GeO2", density=2.196e3, max_temp=1000, abs_coeff=1e-3, n_list_path = ('n_germania.txt',2), k_list_path = ('k_germania.txt',2))
#
# #Print variables.
germania.print_variables()

#Now comment out the code above, uncomment the block below, and re-run the script.
#Here, germania is being loaded from the saved_materials directory.

# germania = Material(name="GeO2")
# germania.print_variables()

#Here's a plot of n and k against wavelength. Uncomment the block below.

# import matplotlib.pyplot as plt
# import numpy as np
# n = []
# k = []
# wavelengths = np.linspace(1e-6, 55e-6, 100)
# for wl in wavelengths:
#     n.append(germania.get_n(wl))
#     k.append(germania.get_k(wl))
# plt.plot(wavelengths*1e6, n, color='b', label='n')
# plt.plot(wavelengths*1e6, k, color='r', label='k')
# plt.xlabel('Wavelength (microns)')
# plt.ylabel('Complex refractive index')
# plt.legend(loc='upper left')
# plt.show()

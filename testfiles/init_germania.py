from Starshot.materials.material import Material

#Initialise germania. Saves GeO2.pkl into saved_materials directory. Note that path is given to a file
#containing refractive index/extinction coefficient vs wavelength.
germania = Material(name="GeO2", density=3.65e3, max_temp=1000, abs_coeff=1e-6, n_list_path = 'n_germania.txt', k_list_path = 'k_germania.txt')

#Print variables.
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
# wavelengths = np.linspace(5, 55, 100)
# for wl in wavelengths:
#     n.append(germania.get_n(wl))
#     k.append(germania.get_k(wl))
# plt.plot(wavelengths, n, color='b', label='n')
# plt.plot(wavelengths, k, color='r', label='k')
# plt.xlabel('Wavelength (microns)')
# plt.ylabel('Complex refractive index')
# plt.legend(loc='upper left')
# plt.show()

from Starshot.materials.material import Material

#Initialise silica. Saves SiO2.pkl into saved_materials directory.
silica = Material(name="SiO2", density=2.196e3, max_temp=1000, abs_coeff=1e-6, n_list_path = None, k_list_path = None)

#Print variables.
silica.print_variables()

#Instead of a table of refractive index vs wavelength, we can add equations. A guide on how to define these equations is written in n_sio2_kkg.py!!!
silica.add_equation(name="n_kkg", range=(1e-6, 25e-6), filepath="n_sio2_kkg.py", n_or_k="n")
silica.add_equation(name="k_kkg", range=(1e-6, 25e-6), filepath="k_sio2_kkg.py", n_or_k="k")

#Note that 'n_kkg' and 'k_kkg' are added to the equations lists.
silica.print_variables()

#Suppose we made a mistake with 'n_kkg' and we wanted to delete it.
silica.rmv_equation(name="n_kkg", n_or_k="n")

#It is now gone.
silica.print_variables()

#Adding it back again...
silica.add_equation(name="n_kkg", range=(1e-6, 25e-6), filepath="n_sio2_kkg.py", n_or_k="n")

#It's back!
silica.print_variables()

#Now comment out the code above, uncomment the code below, and re-run the script.
#Here, silica is being loaded from the saved_materials directory.

# silica = Material(name="SiO2")
# silica.print_variables()

#Here's a plot of n and k against wavelength. Uncomment the block below.

# import matplotlib.pyplot as plt
# import numpy as np
# n = []
# k = []
# wavelengths = np.linspace(1e-6, 25e-6, 100)
# for wl in wavelengths:
#     n.append(silica.get_n(wl))
#     k.append(silica.get_k(wl))
# plt.plot(wavelengths*1e6, n, color='b', label='n')
# plt.plot(wavelengths*1e6, k, color='r', label='k')
# plt.xlabel('Wavelength (microns)')
# plt.ylabel('Complex refractive index')
# plt.legend(loc='upper left')
# plt.show()

import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from sympy.solvers import solve
from sympy import Symbol
from tmm import tmm, general_tmm
from make_transfer_matrix import make_transfer_matrix
from optical_constants import n_silica
from find_eq_temp import find_absorption_from_coefficient, directional_emissivity, spectral_power_flux

# wavelengths = 1e-6*np.array([5,10,15,20,25])
wavelengths = np.linspace(1e-6,25e-6,101)
structure = [(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)]
s2 = [(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)]
s3 = [(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)]
s4 = [(1.45, -206e-9)]
angles = np.linspace(0,pi/2,100)
# ems = []
# plots = []
# for wavelength in wavelengths:
#     ems = []
#     for theta in angles:
#         ems.append(directional_emissivity(theta, structure, wavelength))
#     plots.append(plt.plot(angles*180/(pi), ems)[0])
#
# plt.legend(plots, ("$\lambda$ = 5 $\mu$m","$\lambda$ = 10 $\mu$m","$\lambda$ = 15 $\mu$m","$\lambda$ = 20 $\mu$m","$\lambda$ = 25 $\mu$m"))
# plt.xlim(0,90)
# plt.ylim(0,0.6)
# plt.xlabel("$\\theta$ ($^o$)")
# plt.ylabel("Absorption")
# plt.title("Directional Absorption of A7 structure as a function of angle\nat different wavelengths")
#
# plt.show()

pwr = []
pwr1 = []
pwr2 = []
pwr3 = []
temperature = 300
for wavelength in wavelengths:
    pwr.append(spectral_power_flux(wavelength, s2, temperature))
    pwr1.append(spectral_power_flux(wavelength, structure, temperature))
    pwr2.append(spectral_power_flux(wavelength, s3, temperature))
    pwr3.append(spectral_power_flux(wavelength, s4, temperature))
plt.plot(wavelengths, pwr)
plt.plot(wavelengths, pwr1)
plt.plot(wavelengths, pwr2)
plt.plot(wavelengths, pwr3)

h = 6.62607004e-34       # Planck's constant in SI
c = 299792458             # speed of light in SI
k_B = 1.38064852e-23        # Boltzmann constant in SI
sigma = 5.67e-8
I = ((2*pi*2*h*c**2)/wavelengths**5)*(1/(np.exp(h*c/(wavelengths*k_B*temperature))-1))         # Planck's Law
plt.plot(wavelengths, I)
plt.show()

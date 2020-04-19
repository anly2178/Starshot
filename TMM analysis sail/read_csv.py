# f = open("RefractiveIndexINFO.csv", "r")
# vals = f.readlines()
# for line in vals:
#     if line[-3] == ',':
#         print(line)
#         vals.remove(line)
# f.close()
# g = open("JustK.txt", "w")
# g.writelines(vals)

from optical_constants import n_csi
import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from tmm import tmm, general_tmm
from utility_functions_reflectivity import graph_structure, calc_avg_reflectance
from optical_constants import n_silica
from find_eq_temp import find_absorption_from_coefficient, directional_emissivity, spectral_power_flux
from optical_constants import n_silica, n_csi


# ncsi = (n_csi())
# wl = []
# n = []
# k = []
# for point in ncsi:
#     wl.append(point[0])
#     n.append(point[1].real)
#     k.append(point[1].imag)
#
# wl=np.array(wl)
#
# # plt.plot(wl, n)
# plt.plot(wl, k)
#
# plt.ylabel('Refractive index (n+ik)')
# plt.xlabel('Wavelength ($\mu$m)')
# plt.title('Refractive index of c-Si against wavelength')
# plt.show()

# Check power in at all wavelengths for A5, A7
A1 = [(1.45, -206e-9)]
A3 = [(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)]
A5 = [(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)]
A7 = [(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)]
ratio = 100e9
wavelength_0 = 1.2e-6
betas = (0,0.2)
beta = betas[0]
coeff = np.linspace(1e-6,1e-2,50)
wl = wavelength_0*np.sqrt((1+beta)/(1-beta))

rho_A1 = 0.453
rho_A3 = 0.867
rho_A5 = 1.19
rho_A7 = 1.37

Abs_A5 = []
Abs_A7 = []

    # for wl in wavelength:
A_A1 = find_absorption_from_coefficient(A1, coeff, wl)
A_A3 = find_absorption_from_coefficient(A3, coeff, wl)
A_A5 = find_absorption_from_coefficient(A5, coeff, wl)
A_A7 = find_absorption_from_coefficient(A7, coeff, wl)
P_A1 = A_A1*ratio*rho_A1*(1-beta)/(1+beta)
P_A3 = A_A3*ratio*rho_A3*(1-beta)/(1+beta)
P_A5 = A_A5*ratio*rho_A5*(1-beta)/(1+beta)
P_A7 = A_A7*ratio*rho_A7*(1-beta)/(1+beta)
b1 = plt.plot(coeff, P_A1, 'b')
plt.plot(coeff, P_A3,'b.')
plt.plot(coeff, P_A5,'b-.')
plt.plot(coeff, P_A7,'b--')
print(wl)
beta = betas[1]
wl = wavelength_0*np.sqrt((1+beta)/(1-beta))
print(wl)
A_A1 = find_absorption_from_coefficient(A1, coeff, wl)
A_A3 = find_absorption_from_coefficient(A3, coeff, wl)
A_A5 = find_absorption_from_coefficient(A5, coeff, wl)
A_A7 = find_absorption_from_coefficient(A7, coeff, wl)
P_A1 = A_A1*ratio*rho_A1*(1-beta)/(1+beta)
P_A3 = A_A3*ratio*rho_A3*(1-beta)/(1+beta)
P_A5 = A_A5*ratio*rho_A5*(1-beta)/(1+beta)
P_A7 = A_A7*ratio*rho_A7*(1-beta)/(1+beta)
b2 = plt.plot(coeff, P_A1, 'r')
plt.plot(coeff, P_A3,'r.')
plt.plot(coeff, P_A5,'r-.')
plt.plot(coeff, P_A7,'r--')

# plt.xscale('log')
print(P_A7)

plt.xlabel('absorption coefficient of silica (cm$^{-1}$)')
plt.ylabel('Absorbed power flux from laser (W m$^{-2}$))')
plt.legend((b1[0],b2[0]),('beta = 0', 'beta = 0.2'))
plt.show()

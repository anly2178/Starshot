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
from optical_constants import n_silica

import time
start_time = time.time()


# print(vals)
# print(vals[0]*np.conj(vals[0])+vals[1]*np.conj(vals[1]))
# print(vals[3]*np.conj(vals[3])+vals[2]*np.conj(vals[2]))
# print(((vals[0]*np.conj(vals[0]))+(vals[2]*np.conj(vals[2])))/2)
# print(((vals[1]*np.conj(vals[1]))+(vals[3]*np.conj(vals[3])))/2)
#
# vals2 = tmm([(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)], 1.2e-6)
# print(((vals2[0]*np.conj(vals2[0]))))

angles = np.linspace(0,pi/2,500)
Rp = []
Tp = []
Rs = []
Ts = []
A = []
for theta in angles:
    wl = 100e-9
    k =(n_silica(wl))
    print(k.real)
    # print(k)
    # k = 2e-1
    vals = general_tmm([(k, -100e-9)], wl, theta)
    # vals = general_tmm([(1.45+k*1j, -156e-9), (1, -452e-9), (1.45+k*1j, -161e-9), (1, -448e-9), (1.45+k*1j, -161e-9), (1, -452e-9), (1.45+k*1j, -156e-9)], 1.2e-6, theta)
    # r = (vals[0]*np.conj(vals[0])+vals[2]*np.conj(vals[2]))/2
    # t = (vals[1]*np.conj(vals[1])+vals[3]*np.conj(vals[3]))/2
    rp = vals[0]*np.conj(vals[0])
    # print(rp)
    tp = vals[1]*np.conj(vals[1])
    rs = vals[2]*np.conj(vals[2])
    ts = vals[3]*np.conj(vals[3])

    # print(r)
    Rp.append(rp)
    Tp.append(tp)
    Rs.append(rs)
    Ts.append(ts)
    # print(len(Rp))
    # A.append(1-r-t)
a=plt.plot(angles*180/pi, Rp, 'b')
b=plt.plot(angles*180/pi, Tp, 'b--')
c=plt.plot(angles*180/pi, Rs, 'g')
d=plt.plot(angles*180/pi, Ts, 'g--')
# c=plt.plot(angles, A)
plt.legend((a[0],b[0],c[0],d[0]), ('Rp', 'Tp', 'Rs', 'Ts'))
plt.grid()
plt.xlabel('$\\theta$ ($^o$)')
plt.ylabel('Reflectance/Transmittance')
plt.title('Reflectance and transmittance of s and p-polarised\nincident light, calculated using TMM')
plt.show()

# thickness = np.linspace(100e-9, 1000e-9,100)
# R = []
# for d in thickness:
#     R.append(tmm([(1.45, d)], 400e-9)[0]*np.conj(tmm([(1.45, d)], 400e-9)[0]))
#
# plt.plot(thickness,R)
# plt.show()

# seems like it works!

# n=[]
# k=[]
# wavelength = np.linspace(1e-6, 25e-6, 101)
# for i in wavelength:
#     idx = n_silica(i)
#     n.append(idx.real)
#     k.append(idx.imag)
#
# plt.plot(wavelength*1e6, n)
# plt.plot(wavelength*1e6, k)
#
# plt.ylabel('Refractive index (n+ik)')
# plt.xlabel('Wavelength ($\mu$m)')
# plt.title('Refractive index of SiO$_2$ against wavelength')
# plt.show()

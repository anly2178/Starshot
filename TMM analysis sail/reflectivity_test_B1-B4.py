import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from tmm import tmm
from make_transfer_matrix import make_transfer_matrix
from utility_functions_reflectivity import graph_structure, calc_avg_reflectance, calc_avg_transmittance

""" A recreation of Fig S1, but only B1-B4 """

B2 = graph_structure([(1.45, -5e-9), (3.5, -61e-9)],'b')
B2d = graph_structure([(1.45, -63e-9),(3.5, -54e-9)], 'b--')
B3 = graph_structure([(1.45, -5e-9), (1, -550e-9), (3.5, -61e-9)], 'r')
B3d = graph_structure([(1.45, -73e-9), (1, -619e-9), (3.5, -51e-9)], 'r--')
B4 = graph_structure([(1.45, -5e-9), (3.5, -33e-9), (1, -506e-9), (3.5, -34e-9)], 'g')
B4d = graph_structure([(1.45, -6e-9), (3.5, -31e-9), (1, -523e-9), (3.5, -33e-9)], 'g--')
B4dd = graph_structure([(1.45, -65e-9), (3.5, -45e-9), (1, -631e-9), (3.5, -5e-9)], 'g-.')
B3=graph_structure([(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)])
B4=graph_structure([(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)])

plt.ylim(0,1)
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Reflectivity')
plt.title('Spectral reflectance at normal incidence (recreation of Fig S1)')
plt.legend((B2,B2d,B3,B3d,B4,B4d,B4dd),("B2","B2'","B3","B3'","B4","B4'","B4''"))
plt.show()
#
print('<R> for B2 = '+str(calc_avg_reflectance([(1.45, -5e-9), (3.5, -61e-9)]).real))
print("<R> for B2' = "+str(calc_avg_reflectance([(1.45, -63e-9),(3.5, -54e-9)]).real))
print("<R> for B3 = "+str(calc_avg_reflectance([(1.45, -5e-9), (1, -550e-9), (3.5, -61e-9)]).real))
print("<R> for B3' = "+str(calc_avg_reflectance([(1.45, -73e-9), (1, -619e-9), (3.5, -51e-9)]).real))
print("<R> for B4 = "+str(calc_avg_reflectance([(1.45, -5e-9), (3.5, -33e-9), (1, -506e-9), (3.5, -34e-9)]).real))
print("<R> for B4' = "+str(calc_avg_reflectance([(1.45, -6e-9), (3.5, -31e-9), (1, -523e-9), (3.5, -33e-9)]).real))
print("<R> for B4'' = "+str(calc_avg_reflectance([(1.45, -65e-9), (3.5, -45e-9), (1, -631e-9), (3.5, -5e-9)]).real))

# print(calc_avg_transmittance([(1.45 + 20j, 5e-9), (3.5, 61e-9)]))
# print(calc_avg_transmittance([(1.45, 63e-9), (3.5, 54e-9)]))
# print(calc_avg_transmittance([(1.45, 5e-9), (3.5, 54e-9)]))
#
print(calc_avg_reflectance([(1.45+0.001j, -50e-9)]) + calc_avg_transmittance([(1.45+0.001j, -50e-9)]))

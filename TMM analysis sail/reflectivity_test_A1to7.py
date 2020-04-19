import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from tmm import tmm
from make_transfer_matrix import make_transfer_matrix
from utility_functions_reflectivity import graph_structure, calc_avg_reflectance

""" A recreation of Fig S1, but only A1-A7 """

# Takes the tuple describing the structure and plots it (100 points linspaced)
# over the bandwidth

a=graph_structure([(1.45, -206e-9)])
b=graph_structure([(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)])
c=graph_structure([(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)])
d=graph_structure([(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)])

plt.ylim(0,1)
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Reflectivity')
plt.title('Spectral reflectance at normal incidence (recreation of Fig S1)')
plt.legend((a,b,c,d),('A1','A3','A5','A7'))
plt.show()

print('<R> for A1 = '+str(calc_avg_reflectance([(1.45, -206e-9)]).real))
print('<R> for A3 = '+str(calc_avg_reflectance([(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)]).real))
print('<R> for A5 = '+str(calc_avg_reflectance([(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)]).real))
print('<R> for A7 = '+str(calc_avg_reflectance([(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)]).real))

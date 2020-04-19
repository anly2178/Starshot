import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from find_W import W, rho_S
from find_eq_temp import solveTemp, solveTemp_csi

from matplotlib.pyplot import figure
figure(num=None, figsize=(6.5, 6.5), dpi=80, facecolor='w', edgecolor='k')

abs = np.logspace(-2,-6,15)
plt.xscale('log')

p100_1=solveTemp(100e9,[(1.45, -206e-9)], rho_S(2.196e6,206e-9))
plt.plot(abs, p100_1[0], 'b')
p100_2=solveTemp(100e9,[(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)], rho_S(2.196e6,394e-9))
plt.plot(abs, p100_2[0], 'b')
p100_3=solveTemp(100e9,[(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)], rho_S(2.196e6,542e-9))
plt.plot(abs, p100_3[0], 'b')
p100_4=solveTemp(100e9,[(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)], rho_S(2.196e6,634e-9))
plt.plot(abs, p100_4[0], 'b')

p1000_1=solveTemp(1000e9,[(1.45, -206e-9)], rho_S(2.196e6,206e-9))
plt.plot(abs, p1000_1[0], 'r')
p1000_2=solveTemp(1000e9,[(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)], rho_S(2.196e6,394e-9))
plt.plot(abs, p1000_2[0], 'r')
p1000_3=solveTemp(1000e9,[(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)], rho_S(2.196e6,542e-9))
plt.plot(abs, p1000_3[0], 'r')
p1000_4=solveTemp(1000e9,[(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)], rho_S(2.196e6,634e-9))
plt.plot(abs, p1000_4[0], 'r')

p10_1=solveTemp(10e9,[(1.45, -206e-9)], rho_S(2.196e6,206e-9))
plt.plot(abs, p10_1[0], 'g')
p10_2=solveTemp(10e9,[(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)], rho_S(2.196e6,394e-9))
plt.plot(abs, p10_2[0], 'g')
p10_3=solveTemp(10e9,[(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)], rho_S(2.196e6,542e-9))
plt.plot(abs, p10_3[0], 'g')
p10_4=solveTemp(10e9,[(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)], rho_S(2.196e6,634e-9))
plt.plot(abs, p10_4[0], 'g')

def list_writer(f, list):
    for row in list:
        f.write("[")
        for element in row:
            f.write(str(element))
            f.write(', ')
        f.write(']\n')

f = open('Equilibrium Temperatures', 'w')
f.write("10 GW/g, A1:\n")
list_writer(f, p10_1)
f.write("10 GW/g, A3:\n")
list_writer(f, p10_2)
f.write("10 GW/g, A5:\n")
list_writer(f, p10_3)
f.write("10 GW/g, A7:\n")
list_writer(f, p10_4)
f.write("\n==========================================\n")
f.write("100 GW/g, A1:\n")
list_writer(f, p100_1)
f.write("100 GW/g, A3:\n")
list_writer(f, p100_2)
f.write("100 GW/g, A5:\n")
list_writer(f, p100_3)
f.write("100 GW/g, A7:\n")
list_writer(f, p100_4)
f.write("\n==========================================\n")
f.write("1000 GW/g, A1:\n")
list_writer(f, p1000_1)
f.write("1000 GW/g, A3:\n")
list_writer(f, p1000_2)
f.write("1000 GW/g, A5:\n")
list_writer(f, p1000_3)
f.write("1000 GW/g, A7:\n")
list_writer(f, p1000_4)
f.close()


plt.ylim(0,1850)
plt.xlim(10**(-6),10**(-2))
plt.yticks((0,500,1000,1500))
plt.grid(True)

plt.title('Equibrium temperatures for structures A1-A7 for different\npower-mass ratios against absorption coefficient (Fig 2d)')
plt.xlabel('Absorption coefficient (cm$^{-1}$)')
plt.ylabel('T (K)')
plt.show()

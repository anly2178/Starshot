import matplotlib.pyplot as plt
import numpy as np
#
# absorption = np.logspace(-9, -4, 100) #absolute absorption
# I = 1e10 #Wm^-2, the intensity of the laser
# P_in = absorption * I * 10 #W, Power absorbeds
# sigma = 5.6704e-8 #Wm^-2K^-4, Stefan-Boltzmann constant
# T_space = 2.73 #K, temperature in space
# emissivity = [0.5, 0.1, 1e-2, 1e-3] #Emissivity
#
# i = 0
# linestyles = ['-','--','-.', ':']
# while i < len(emissivity):
#     T = (P_in / (20*emissivity[i]*sigma) + T_space**4)**0.25 #using Stefan-Boltzmann Law
#     plt.plot(absorption, T, linestyles[i], label = 'emissivity = {}'.format(emissivity[i]))
#     i += 1
#
#
# plt.suptitle('Fig 4b')
# plt.ylim(200,2000)
# plt.ylabel('Equilibrium Temperature (K)')
# plt.xlim(1e-9, 1e-4)
# plt.xlabel('Absorption')
# plt.xscale('log')
# plt.legend(loc = 'lower right')
# plt.show()

absorb_coeff = np.logspace(-4,0,100) #m-1
emissivity = 0.01 #absolute emissivity
sigma = 5.6704e-8 #Wm^-2K^-4, Stefan-Boltzmann constant
P0_on_mp = [1e13, 1e14, 1e15] #W/kg
thickness_A1 = 206e-9 #m
rho_A1 = 0.453e-3 #kg/m^2
i = 0

while i < 3:
    P0_mp = P0_on_mp[i]
    T = (P0_mp * rho_A1 * (1-np.exp(-absorb_coeff * thickness_A1)) / (emissivity * sigma))**0.25
    plt.plot(absorb_coeff, T)
    i += 1

plt.ylim(0,1800)
plt.ylabel('Equilibrium Temperature (K)')
plt.xlim(1e-4, 1)
plt.xlabel('Absorption coefficient [m^-1]')
plt.xscale('log')
plt.show()

"""
This is a template showing how to use the scripts.
NOTE: This file would be run in the parent directory.
"""
#Import modules
from starshot_lib.temp import *
from starshot_lib.motion import *
import numpy as np
import matplotlib.pyplot as plt
from starshot_lib.write_data import write_data

#If absorption coefficient is known and the material is available, one can calculate
#the absorptance using the function, otherwise it must be calculated manually and entered.
abs_coeff = 1e-3
material = 'SiO2'
thickness = 1e-6
wavelength = 1064e-9
A = find_absorptance(material, thickness, abs_coeff, wavelength)

#Initialise parameters of sail and DE system. These terms are defined in README.md
params = {"m_sail": 1e-3, "thickness": 1e-6, "density": 1400, "reflectivity": 1, "absorptance": A,
            "k": 1, "power": 1e11, "laser_size": 1e4, "wavelength": 1064e-9, "alpha": 1}

#Calculate the state over time.
state, time = state_vs_t(params)

#Stores results in tabulated form in txt file
filepath = 'sample.txt'
write_data(params, state, time, filepath)

#Extract rows of state, which include speed (as fraction of c), distance, temperature
betas = state[0,:]
dist = state[1,:]
temp = state[2,:]

#Plots temperature vs time and saves it
plt.figure(1)
plt.plot(t, T1)
plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Equilibrium temperature (K)')
plt.legend()
plt.savefig('sample_tempvst.png')
#Plots temperature vs beta and saves it
plt.figure(2)
plt.plot(betas1, T1)
plt.xlabel('Beta')
plt.ylabel('Equilibrium temperature (K)')
plt.legend()
plt.savefig('sample_tempvsbeta.png')
#Plots temperature vs distance and saves it
plt.figure(3)
plt.plot(dist1, T1)
plt.xscale('log')
plt.xlabel('Distance (m)')
plt.ylabel('Equilibrium temperature (K)')
plt.legend()
plt.savefig('sample_tempvsdist.png')

#Shows the figures.
plt.show()

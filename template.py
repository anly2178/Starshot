"""
This is a template showing how to use the scripts.
NOTE: This file would be run in the parent directory.

One could easily move this file into the parent directory,
replace the values in here and use this as a 'calculator'.
"""
#Import modules
from starshot_lib.motion import *
from starshot_lib.temp import add_all_temp
from starshot_lib.fill_params import fill_params
import numpy as np
import matplotlib.pyplot as plt
from starshot_lib.write_data import write_data

#Initialise parameters of sail and DE system. These terms are defined in README.md
params = {  "material": 'SiO2',
            "m_sail": 1e-3,
            "thickness": None,
            "area": 10,
            "density": 1400,
            "abs_coeff": 1e-3,
            "absorptance": None,
            "reflectance": None,
            "transmittance": None,
            "power": 1e11,
            "diameter": 1e4,
            "wavelength": 1064e-9}

#Fills the missing parameters, if possible - in this case, thickness and absorptance
fill_params(params)

#Calculate the state over time.
state, time = state_vs_t(params)

#Includes the temperature at each state
state = add_all_temp(state)

#Stores results in tabulated form in txt file
filepath = 'sample.txt'
write_data(params, state, time, filepath)

#Extract rows of state, which include speed (as fraction of c), distance, temperature
betas = state[0,:]
dist = state[1,:]
temp = state[2,:]

#Plots temperature vs time and saves it
plt.figure(1)
plt.plot(time, temp)
plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Equilibrium temperature (K)')
plt.savefig('sample_tempvst.png')
#Plots temperature vs beta and saves it
plt.figure(2)
plt.plot(betas, temp)
plt.xlabel('Beta')
plt.ylabel('Equilibrium temperature (K)')
plt.savefig('sample_tempvsbeta.png')
#Plots temperature vs distance and saves it
plt.figure(3)
plt.plot(dist, temp)
plt.xscale('log')
plt.xlabel('Distance (m)')
plt.ylabel('Equilibrium temperature (K)')
plt.savefig('sample_tempvsdist.png')

#Shows the figures.
plt.show()

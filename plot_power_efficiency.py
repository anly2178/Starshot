#Plots Kulkarni 2018 figure 3.
#Shows the efficiency of laser power transfer as a function of velocity.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

#Prepare variables
betas = np.linspace(0,1,100)
pi_forward = betas
pi_sail = np.zeros(100)
pi_backward = np.zeros(100)

#Solve for the variables as a function of beta
i = 0
while i < 100:
    beta = betas[i]
    pi_sail[i] = 2*beta*(1-beta)/(1+beta)
    pi_backward[i] = (1-beta)**2 / (1+beta)
    i += 1

#Plot results
fig, ax = plt.subplots()
ax.plot(betas, pi_forward,':', label = r'$\pi_{\rightarrow}$')
ax.plot(betas, pi_backward,'--', label = r'$\pi_{\leftarrow}$')
ax.plot(betas, pi_sail, label = r'$\pi_{sc}$')

#Plot design
ax.tick_params(which = 'both', direction = 'in', left = True, right = True, bottom = True, top = True)
ax.tick_params(labelleft = True, labelright = False, labelbottom = True, labeltop = False)
ax.set_xlim(0,1)
ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.set_ylim(0,1)
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.set_xlabel(r'$\beta$', fontsize = 13)
ax.legend()

plt.show()

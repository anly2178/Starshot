import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from tmm import tmm
from make_transfer_matrix import make_transfer_matrix

def dW(beta, structure, rho_S):
    gamma = 1/np.sqrt(1-beta**2)
    wavelength = 1.2e-6*np.sqrt((1+beta)/(1-beta))
    dW = np.sqrt(rho_S)/(tmm(structure, wavelength)[0]*np.conj(tmm(structure, wavelength)[0])).real * (gamma*beta)/(1-beta)**2
    return dW

# Function for evaulating RAAD as the integral of dW
def W(structure, rho_S):
    W = integrate.quad(dW, 0, 0.2, args=(structure, rho_S))
    return W

# Function for finding area density from mass density and layer width (nm)
def rho_S(rho, width):
    return rho*width

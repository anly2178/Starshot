import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from tmm import tmm
from make_transfer_matrix import make_transfer_matrix

def graph_structure(structure, color = ''):
    bandwidth = np.linspace(1.2, 1.225*1.2, 100)
    i = 0
    yvals = []
    while i < len(bandwidth):
        r = tmm(structure, bandwidth[i]*1e-6)[0]
        R = (r*np.conj(r)).real
        yvals.append(R)
        i += 1
    p = plt.plot(bandwidth, yvals, color)
    return p[0]

def calc_avg_reflectance(structure):
    sum = 0
    bandwidth = np.linspace(1.2, 1.225*1.2, 100)
    i = 0
    while i < len(bandwidth):
        r = tmm(structure, bandwidth[i]*1e-6)[0]
        R = r*np.conj(r)
        sum += R
        i += 1
    return sum/100

def calc_avg_transmittance(structure):
    sum = 0
    bandwidth = np.linspace(1.2, 1.225*1.2, 100)
    i = 0
    while i < len(bandwidth):
        r = tmm(structure, bandwidth[i]*1e-6)[1]
        R = r*np.conj(r)
        sum += R
        i += 1
    return sum/100

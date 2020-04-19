import scipy
import numpy as np
from numpy import sin, cos, pi

""" Matrix params are a single tuple of size 2. Format (n, d) where n is
    the refractive index, and d is the thickness of the layer. Returns the
    correct transfer matrix. Also needs a wavenumber as input (calculated from
    tmm to reduce computational time).
"""

def make_transfer_matrix(matrix_params, wavenumber):
    n = matrix_params[0]
    d = matrix_params[1]
# Remark that assuming there is no phase change between each layer, the
# wavenumbers at the mth layer is just n_m*k0
    k = wavenumber*n
    M_11 = cos(k*d)
    M_12 = sin(k*d)/k
    M_21 = -k*sin(k*d)
    M_22 = cos(k*d)
    M = [[M_11, M_12], [M_21, M_22]]
    return M

def make_p_transfer_matrix(matrix_params, wavenumber, angle):
    n = matrix_params[0]
    d = matrix_params[1]
    k = wavenumber*n

# Other values in terms of previous parameters
    delta = k*d*cos(angle)
    eng = 1j*k/cos(angle)       # Note: not really eng in McLeod, but should be the analogue for this analysis. To get theirs, just divide by impedance of free space

    M_11 = cos(delta)
    M_12 = 1j*sin(delta)/eng
    M_21 = 1j*eng*sin(delta)
    M_22 = cos(delta)
    M = [[M_11, M_12], [M_21, M_22]]
    return M

def make_s_transfer_matrix(matrix_params, wavenumber, angle):
    n = matrix_params[0]
    d = matrix_params[1]
    k = wavenumber*n

# Other values in terms of previous parameters
    delta = k*d*cos(angle)
    eng = 1j*k*cos(angle)       # Note: not really eng in McLeod, but should be the analogue for this analysis. To get theirs, just divide by impedance of free space

    M_11 = cos(delta)
    M_12 = 1j*sin(delta)/eng
    M_21 = 1j*eng*sin(delta)
    M_22 = cos(delta)
    M = [[M_11, M_12], [M_21, M_22]]
    return M

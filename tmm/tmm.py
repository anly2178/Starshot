import scipy
import numpy as np
from numpy import sin, cos, pi, arcsin
from make_transfer_matrix import make_transfer_matrix, make_p_transfer_matrix, make_s_transfer_matrix

""" This function will take a set of transfer matrices and return in a double
    (r,t), the reflectivity and transmittance coefficients. From these,
    need to take the modulus squared to find Reflectivity and Transmittance.
    Arguments are the parameters based on the layered structures in a
    tuple, as well as the wavelength of light this is being calculated for.

    In more detail, matrix_params = ( (n_1, d_1), (n_2, d_2), ... , (n_m, d_m) ),
    A TUPLE OF TUPLES
    With the mth layer being the closest to the beam, and 1st layer the furthest.

    MAKE SURE TO ENTER NEGATIVE DISTANCES IN MATRICES - this will affect r,t
    especially for complex refractive indices

    Take theta to be the initial angle of incidence, which will
    change as you go down layers due to refraction. Remember theta in radians.
    Notably, we need to consider p and s polarisations, so we have two r's and
    two t's which will be returned
"""

def tmm(matrix_params, wavelength, theta):
    n0 = 1                  # refractive index of vacuum
    k0 = n0*2*pi/wavelength     # WAVELENGTH IN METRES
    p_field_vector = [[1],[1j*k0/cos(theta)]]
    s_field_vector = [[1],[1j*k0*cos(theta)]]
    M_p = np.array([[1,0],[0,1]])
    M_s = np.array([[1,0],[0,1]])

    # Here, we need to just make the matrices and continue to left-multiply them
    # to the field vector, beginning with M_1, ending at M_m
    for matrix_param in matrix_params:
        new_theta = arcsin(n0/matrix_param[0]*sin(theta))       # Simple and speedy
        M_p = np.matmul(make_p_transfer_matrix(matrix_param, k0, new_theta), M_p)
        M_s = np.matmul(make_s_transfer_matrix(matrix_param, k0, new_theta), M_s)

    p_field_vector = np.matmul(M_p,p_field_vector)
    s_field_vector = np.matmul(M_s,s_field_vector)
    # The field vector is now the col vector (alpha,beta), which we can use to
    # find r, t:

    # Note we can use theta here since angle of incidence = angle of outgoing ray
    E_p = p_field_vector[0, 0]
    H_p = p_field_vector[1, 0]
    H_p_divided = H_p/(1j*k0/cos(theta))     # For easier calcs

    E_s = s_field_vector[0, 0]
    H_s = s_field_vector[1, 0]
    H_s_divided = H_s/(1j*k0*cos(theta))     # For easier calcs

# Calc coefficients for p-polarised light
    a_i_p = (E_p+H_p_divided)/2
    a_r_p = (E_p-H_p_divided)/2

    r_p = a_r_p/a_i_p
    t_p = 1/a_i_p

# Calc coefficients for s-polarised light
    a_i_s = (E_s+H_s_divided)/2
    a_r_s = (E_s-H_s_divided)/2

    r_s = a_r_s/a_i_s
    t_s = 1/a_i_s

    return (r_p,t_p,r_s,t_s)

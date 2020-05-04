import numpy as np
from .temp import find_temp
from .laser import find_fraction_incident

def differential_eq(x, params):
    """
    Defines the differential equations as in Kulkarni 2018 equation (23).
    Relates beta and distance to time t in seconds.
    Assume Gaussian beam with diffraction effects.
    Assume optimal mass condition: sail_mass = payload_mass.
    x is a state vector containing beta and the position/distance.
    """
    #Get state of sail at that time
    beta = x[0]
    dist = x[1]
    #Get parameters
    m_sail = params["m_sail"]
    m_tot = 2*m_sail
    reflectance = params["reflectance"]
    power = params["power"]

    fraction = find_fraction_incident(params, dist) #fraction of power incident
    lor = 1/(1-beta**2)**0.5 #Lorentz factor
    #Derivative of state with respect to time
    beta_dot = 2 * reflectance * power * fraction * (1-beta) / (m_tot * c**2 * lor**3 * (1+beta))
    vel = c*beta
    xdot = np.array([beta_dot, vel])
    return xdot

def state_vs_t(params):
    """
    Returns an array of states
        First row contains the speed of sail as a fraction of the speed of light
        Second row contains the distance from the DE system
    and an array of corresponding times.
    Assumes diffraction effects.

    Parameters must be defined (in SI units) and passed into function. For example:
        params = {"m_sail": 1e-3, "thickness": 1e-6, "density": 1400, "reflectance": 1, "absorptance": 9e-8,
                    "power": 1e11, "diameter": 1e4, "wavelength": 1064e-9}
    """
    c = 2.998e8

    #Initialise conditions
    f = lambda t, x : differential_eq(x, params) #Differential equation
    x0 = np.array([0,0])  #Initial state
    t0 = 0  #Initial t
    tf = 1e4  #Final t
    dt = 1  #Step size

    #Create arrays to fill
    t = np.append(np.linspace(0,0.8,5), np.logspace(0,4,195)) #Create time, starts off linear and transitions into logarithmic
    nt = t.size
    nx = x0.size
    x = np.zeros((nx,200))
    x[:,0] = x0 #Speed and distance

    #Runge Kutta method to find states as a function of time
    for i in range(nt - 1):
        dt = t[i+1] - t[i]

        #Find next beta
        k1 = dt * f(t[i], x[:,i])
        k2 = dt * f(t[i] + dt/2, x[:,i] + k1/2)
        k3 = dt * f(t[i] + dt/2, x[:,i] + k2/2)
        k4 = dt * f(t[i] + dt, x[:,i] + k3)
        dx = (k1 + 2*k2 + 2*k3 + k4)/6
        x[:,i+1] = x[:,i] + dx

    return x, t

def grating_differential_eq(x, params, coeffs, angles):
    """
    Defines the differential equations as in Kulkarni 2018 equation (23).
    Changed to account for reflection (and transmission) at an angle.
    Assume Gaussian beam with diffraction effects.
    Relates beta and distance to time t in seconds.
    Assume optimal mass condition: sail_mass = payload_mass.
    Input:  x -- a state vector containing beta and the position/distance.
            params -- dictionary of parameters defining the sail/laser
            coeffs -- an array of reflection and transmission coefficients
            angles -- array of angles corresponding to the reflected/transmitted orders
    """
    c = 2.998e8 #ms^-1
    #Get state of sail at that time
    beta = x[0]
    dist = x[1]
    #Get parameters
    m_sail = params["m_sail"]
    m_tot = 2*m_sail
    reflectance = params["reflectance"]
    power = params["power"]

    fraction = find_fraction_incident(params, dist) #fraction of power incident
    lor = 1/(1-beta**2)**0.5 #Lorentz factor
    #Derivative of state with respect to time
    beta_dot = 2 * reflectance * power * fraction * (1-beta) / (m_tot * c**2 * lor**3 * (1+beta))
    vel = c*beta
    xdot = np.array([beta_dot, vel])
    return xdot

def grating_state_vs_t(params):
    """
    For grating, returns an array of states
        First row contains the speed of sail as a fraction of the speed of light
        Second row contains the distance from the DE system
    and an array of corresponding times.
    Assumes diffraction effects.

    Parameters must be defined (in SI units) and passed into function. For example:
        params = {"m_sail": 1e-3, "thickness": 1e-6, "density": 1400, "reflectance": 1, "absorptance": 9e-8,
                    "power": 1e11, "diameter": 1e4, "wavelength": 1064e-9}
    """

    #Initialise conditions
    f = lambda t, x : differential_eq(x, params) #Differential equation
    x0 = np.array([0,0])  #Initial state
    t0 = 0  #Initial t
    tf = 1e4  #Final t
    dt = 1  #Step size

    #Create arrays to fill
    t = np.append(np.linspace(0,0.8,5), np.logspace(0,4,195)) #Create time, starts off linear and transitions into logarithmic
    nt = t.size
    nx = x0.size
    x = np.zeros((nx,200))
    x[:,0] = x0 #Speed and distance

    #Runge Kutta method to find states as a function of time
    for i in range(nt - 1):
        dt = t[i+1] - t[i]

        #Find next beta
        k1 = dt * f(t[i], x[:,i])
        k2 = dt * f(t[i] + dt/2, x[:,i] + k1/2)
        k3 = dt * f(t[i] + dt/2, x[:,i] + k2/2)
        k4 = dt * f(t[i] + dt, x[:,i] + k3)
        dx = (k1 + 2*k2 + 2*k3 + k4)/6
        x[:,i+1] = x[:,i] + dx

    return x, t

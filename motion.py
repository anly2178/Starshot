import numpy as np

def differential_eq(x, sail):
    """Returns acceleration and speed of lightsail.

    Parameters
    ----------
    list of floats
        x is a list containing beta (v/c) in the 0th index and distance [m] in
        the 1th index. i.e. [beta, dist]
    Sail
        Instance of Sail class

    Returns
    -------
    list of floats
        Returns a list containing the acceleration and speed. [beta_dot, speed [m/s]]
    """
    c = 2.998e8 #ms^-1
    #Get state of sail at that time
    beta, dist = x
    #Get parameters
    tot_mass = 2 * sail.mass #Optimal mass condition
    angles_coeffs = sail.angles_coeffs #degrees
    fraction = sail._find_fraction(dist) #fraction of power incident
    power_inc = fraction * sail.power #W

    power_ref = 0
    for angle, r, t in angles_coeffs:
        fac = np.cos(np.deg2rad(angle))
        p_r = r * power_inc * fac
        p_t = t * power_inc * (1 - fac)
        power_ref += p_r + p_t

    lor = 1/(1-beta**2)**0.5 #Lorentz factor
    #Derivative of state with respect to time
    beta_dot = 2 * power_ref * (1-beta) / (tot_mass * c**2 * lor**3 * (1+beta))
    vel = c*beta
    xdot = np.array([beta_dot, vel])
    return xdot

def state_vs_t(sail):
    """Returns speed/distance array and corresponding time array.

    Parameters
    ----------
    Sail
        Instance of Sail class

    Returns
    -------
    array of arrays of floats
        [[beta, distance [m]], ...]
    array of floats
        [time [s], ...]
    """
    #Initialise conditions
    f = lambda t, x : differential_eq(x, sail) #Differential equation
    x0 = np.array([0,0])  #Initial state
    t = np.append(np.linspace(0,0.8,10), np.logspace(0,4,140)) #Create time, starts off linear and transitions into logarithmic
    x = np.zeros((x0.size,t.size))
    x[:,0] = x0 #Speed and distance
    #Runge Kutta method to find states as a function of time
    for i in range(t.size - 1):
        dt = t[i+1] - t[i]
        #Find next beta
        k1 = dt * f(t[i], x[:,i])
        k2 = dt * f(t[i] + dt/2, x[:,i] + k1/2)
        k3 = dt * f(t[i] + dt/2, x[:,i] + k2/2)
        k4 = dt * f(t[i] + dt, x[:,i] + k3)
        dx = (k1 + 2*k2 + 2*k3 + k4)/6
        x[:,i+1] = x[:,i] + dx
    return x, t

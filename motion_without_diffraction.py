import numpy as np

def no_diff_beta_dot(beta):
    """
    Defines the differential equation in Kulkarni 2018, equation (10).
    Relates beta (= fraction of speed of light) to t/t_rel, where
    t_rel = mc^2 / P.
    Assumes no diffraction effects.
    """
    beta_dot = 2*(1-beta[0])*(1 - beta[0]**2)**1.5 /(1+beta[0])
    return beta_dot

def no_diff_state_vs_t():
    """
    Returns an array of states
        First row contains the speed of sail as a fraction of the speed of light
        Second row contains the distance from the DE system
    and an array of corresponding times.
    Assumes all photons reflect off the sail.
    """
    c = 2.998e8

    #Initialise conditions
    f = lambda t, beta : no_diff_beta_dot(beta) #Differential equation
    x0 = np.array([0,0])  #Initial state
    t0 = 0  #Initial t/t_rel
    tf = 2.5  #Final t/t_rel
    dt = 1e-4  #Initial step size

    #Create arrays to fill
    t = np.arange(t0, tf, dt)
    nt = t.size
    nx = x0.size
    x = np.zeros((nx,nt))
    x[:,0] = x0

    #Runge Kutta method to find states as a function of time
    for i in range(nt - 1):
        k1 = dt * f(t[i], x[:,i])
        k2 = dt * f(t[i] + dt/2, x[:,i] + np.array([k1/2,0]))
        k3 = dt * f(t[i] + dt/2, x[:,i] + np.array([k2/2,0]))
        k4 = dt * f(t[i] + dt, x[:,i] + np.array([k3,0]))
        dbeta = (k1 + 2*k2 + 2*k3 + k4)/6

        x[0,i+1] = x[0,i] + dbeta

        velocity = c * (x[0,i+1] + x[0,i])/2 #Trapezoidal rule
        dposition = velocity * dt
        x[1,i+1] = x[1,i] + dposition

    return x, t

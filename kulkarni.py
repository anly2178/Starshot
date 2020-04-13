#Used https://www.youtube.com/watch?v=1FYrnwqWQNY
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def no_diffraction_beta_dot(beta):
    """
    Defines the differential equation as in Kulkarni 2018 equation (10).
    Relates beta (= fraction of speed of light) to t/t_rel, where
    t_rel = mc^2 / P.
    Assumes no diffraction effects
    """

    beta_dot = 2*(1-beta[0])*(1 - beta[0]**2)**1.5 /(1+beta[0])
    return beta_dot

def diffraction_beta_dot(x, params):
    """
    Defines the differential equations as in Kulkarni 2018 equation (23).
    Relates beta (= fraction of speed of light) to time t in seconds.
    Assume diffraction effects past a critical distance.
    Assume optimal mass condition: sail_mass = payload_mass.
    x is a state vector containing beta and the position/distance.
    k is a constant related to the shape of the sail.
        k = 1 for square sail
        k = pi / 4 for circular sail
    alpha is a constant related to the shape of the DE system.
    """
    c = 2.998e8 #ms^-1

    beta = x[0]
    dist = x[1]

    m_sail = params["m_sail"]
    thickness = params["thickness"]
    density = params["density"]
    k = params["k"]
    power = params["power"]
    laser_size = params["laser_size"]
    wavelength = params["wavelength"]
    alpha = params["alpha"]

    sail_size = np.sqrt(m_sail / (k*density*thickness))
    critical_dist = laser_size * sail_size / (2*wavelength*alpha)
    m_tot = 2*m_sail

    if dist <= critical_dist:
        beta_dot = 2 * power * (1-beta**2)**1.5 * (1-beta) / (m_tot * c**2 * (1+beta))
    else:
        beta_dot = 2 * power * (1-beta**2)**1.5 * (1-beta) * critical_dist**2\
        / (m_tot * c**2 * (1+beta) * dist**2)

    return beta_dot

def runge_kutta4(f, x0, t0, tf, dt):
    """
    Implementation of fourth order Runge Kutta method.
    f is the function, x0 is the initial state, t0 is the initial time, dt is the step size.
    """

    t = np.arange(t0, tf, dt)
    nt = t.size
    nx = x0.size
    x = np.zeros((nx,nt))
    x[:,0] = x0

    for i in range(nt - 1):
        k1 = dt * f(t[i], x[:,i])
        k2 = dt * f(t[i] + dt/2, x[:,i] + k1/2)
        k3 = dt * f(t[i] + dt/2, x[:,i] + k2/2)
        k4 = dt * f(t[i] + dt, x[:,i] + k3)

        dx = (k1 + 2*k2 + 2*k3 + k4)/6
        x[:,i+1] = x[:,i] + dx
    return x, t

def plot_beta_vs_t_on_trel():
    """
    Plots the graph for beta as a function of t/t_rel, where
    t_rel = mc^2 / P.
    """
    f = lambda t, beta : no_diffraction_beta_dot(beta) #Differential equation
    beta0 = np.array([0])  #Initial state
    t0 = 0  #Initial t/t_rel
    tf = 2.5  #Final t/t_rel
    dt = tf * 1e-4  #Step size

    #Solution
    beta, t = runge_kutta4(f, beta0, t0, tf, dt)
    print(t)
    print(beta[0,:])

    #Plot relativistic results
    fig, ax = plt.subplots()
    ax.plot(t, beta[0,:], label = 'Relativistic')

    #Plot non-relativistic results
    nonrel_beta = 2 * t
    ax.plot(t, nonrel_beta, '--', label = 'Non-relativistic')

    #Plot design
    ax.tick_params(which = 'both', direction = 'in', left = True, right = True, bottom = True, top = True)
    ax.tick_params(labelleft = True, labelright = False, labelbottom = True, labeltop = False)

    ax.set_xlim(0,2.5)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.set_ylim(0,1)
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    ax.set_ylabel(r'$\frac{v}{c}$', fontsize = 16, rotation = 0)
    ax.set_xlabel('t/t$_{rel}$', fontsize = 13)

    ax.legend()

    plt.show()

def plot_power_changes():
    """
    Plots figure 3 of Kulkarni 2018, which is the graph of the efficiency
    of the transfer of laser power as a function of beta.
    """
    betas = np.linspace(0,1,100)
    pi_forward = betas
    pi_sail = np.zeros(100)
    pi_backward = np.zeros(100)

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

def calculate_state_vs_t(params):
    """
    Implementation of fourth order Runge Kutta method.
    f is the function, x0 is the initial state, t0 is the initial time, dt is the step size.
    Includes adaptive step size for efficiency.
    """
    f = lambda t, x : diffraction_beta_dot(x, params) #Differential equation
    x0 = np.array([0,0])  #Initial state
    t0 = 0  #Initial t/t_rel
    tf = 1e4  #Final t/t_rel
    dt = 1  #Initial step size

    c = 2.998e8

    dt_min = 0.01 #Minimal step size

    #Set change tolerances
    dbeta_max = 0.001
    dbeta_min = 0.0005

    # t = np.arange(t0, tf, dt)
    t = np.logspace(0,4,1000)
    nt = t.size

    nx = x0.size
    x = np.zeros((nx,1000))
    x[:,0] = x0

    for i in range(nt - 1):
        dt = t[i+1] - t[i]

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

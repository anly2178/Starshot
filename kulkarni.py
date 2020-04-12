#Used https://www.youtube.com/watch?v=1FYrnwqWQNY
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def beta_vs_t_on_trel(beta):
    """
    Defines the differential equation as in Kulkarni 2018 equation (10).
    Relates beta (= fraction of speed of light) to t/t_rel, where
    t_rel = mc^2 / P.
    """

    beta_dot = 2*(1-beta)*(1 - beta**2)**1.5 /(1+beta)
    return beta_dot

def runge_kutta4(f, x0, t0, tf, dt):
    """
    Implementation of fourth order Runge Kutta method.
    f is the function, x0 is the initial state, t0 is the initial time, dt is the step size.
    """

    t = np.arange(t0, tf, dt)
    nt = t.size
    x = np.zeros(nt)
    x[0] = x0

    for i in range(nt - 1):
        k1 = dt * f(t[i], x[i])
        k2 = dt * f(t[i] + dt/2, x[i] + k1/2)
        k3 = dt * f(t[i] + dt/2, x[i] + k1/2)
        k4 = dt * f(t[i] + dt, x[i] + k3)

        dx = (k1 + 2*k2 + 2*k3 + k4)/6
        x[i+1] = x[i] + dx
    return x, t

def plot_beta_vs_t_on_trel():
    """
    Plots the graph for beta as a function of t/t_rel, where
    t_rel = mc^2 / P.
    """
    f = lambda t, beta : beta_vs_t_on_trel(beta) #Differential equation
    beta0 = 0  #Initial state
    t0 = 0  #Initial t/t_rel
    tf = 2.5  #Final t/t_rel
    dt = tf * 1e-4  #Step size

    #Solution
    beta, t = runge_kutta4(f, beta0, t0, tf, dt)

    #Plot relativistic results
    fig, ax = plt.subplots()
    ax.plot(t, beta, label = 'Relativistic')

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

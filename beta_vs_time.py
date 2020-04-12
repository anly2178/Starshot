#Used https://www.youtube.com/watch?v=1FYrnwqWQNY
import numpy as np

def beta_vs_time(beta):
    """
    Defines the differential equation as in Kulkarni 2018 equation (10).
    Relates beta (= fraction of speed of light) to t/t_rel, where
    t_rel = mc^2 / P.
    """

    beta_dot = 2*(1-beta)*(1 - beta**2)**1.5 /(1+beta)
    return beta_dot

#Fourth order Runge Kutta method.
def RungeKutta4(f, x0, t0, tf, dt):
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

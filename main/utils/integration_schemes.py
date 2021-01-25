import numpy as np
from scipy.integrate import solve_ivp

def scipy_integrate_solve_ivp(Y, t, dt, f, f_args):
    # only one time step.
    f_ = lambda t, y : f(y, t, *f_args)
    Y_pot = solve_ivp(f_,(t, t+dt), np.array(Y)).y
    sol = Y_pot[:, 1]
    return sol # return the last y computed

def rk4(Y, t, dt, f, f_args):
    k1=f(Y,t,*f_args)
    k2=f(Y+.5*dt*k1, t+.5*dt,*f_args)
    k3=f(Y+.5*dt*k2, t+.5*dt,*f_args)
    k4=f(Y+dt*k3, t+dt,*f_args)
    Y_pot=Y+dt*(1/6*k1+1/3*k2+1/3*k3+1/6*k4)
    return Y_pot

def euler_explicit(Y, t, dt, f, f_args):
    Y_pot = Y+dt*f(Y,t,*f_args)
    return Y_pot

    

from __future__ import print_function, division
import numpy as np
import scipy as sp
from scipy import constants as const
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def M0_dt(X, P, t):
    # TODO map X, P to locals
    # phi_0 == tan**-1(V_ex / V_ey)
    
    xbar_e = x_e + p_gc(cos(phi_0 + omega_0 * dt) - cos(omega_0)
    ybar_e = y_e - p_gc(sin(phi_0 + omega_e * dt) - sin(omega_0)
    
    xbar_i = x_i + (p_ix / m_i) * dt
    ybar_i = y_i + (p_iy / m_i) * dt
    zbar_i = z_i + (p_iz / m_i) * dt
    zbar_e = z_e + (p_ez / m_e) * dt
    
    return (X, P)


def Mc_dt(X, P, t):
    # TODO map X, P to locals
    b = np.sqrt((x_i - x_e)**2 + (y_i - y_e)**2 + (z_i - z_e)**2)
    alpha = (const.e * Q_i / 4 * const.pi * const.epsilon_0) * dt
    
    pbar_ix = p_ix - alpha(x_i - x_e) / b**3
    pbar_iy = p_iy - alpha(y_i - y_e) / b**3
    pbar_iz = p_iz - alpha(z_i - z_e) / b**3
    
    return (X, P)


def M_dt(X, P, t):
    # not sure how to apply as symplectic


def run():
    # calculate the change in the ion momentum
    ion_charge = 1 # proton, may be neg.
    mass_i = constants.m_p
    mass_e = constants.m_e
    
    ts = np.linspace(0, 12, 100)
    
    M_dt = odeint(M_dt, P0, ts)
    
    

if __name__ == '__main__':
    
    magnetic_field = 1. # Tesla, will be user defined
    
    # Gyrotron Frequency
    omega_e = np.abs(const.e * magnetic_field) / const.m_e
    
    dt_max = 1/8 * 2*const.pi / omega_e
    
    N = [1,2,4,8]
    
    dt_vec = [1/_N * dt_max for _N in N]
    

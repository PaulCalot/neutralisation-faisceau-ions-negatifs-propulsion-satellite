# to compare some results with :
# http://fisica.ufpr.br/sharipov/tube.html

from main import get_maxwellian_mean_speed_from_temperature
from main import available_particles
import numpy as np

R = 8.314 # J/mol.K

def speed_of_sound(gamma,T,M):
    """Return the speed of sound.

    Args:
        gamma (float): adiabatic constant 
        T (float): absolute temperature (K)
        M (float): molecular weight (kg/mol)

    Returns:
        float: c - speed of sound
    """
    return np.sqrt(gamma*R*T/M)

T = 300

# Xe
gamma_Xe = 1.65
M_Xe = 0.131 # kg/mol
m_Xe = available_particles['Xe']['mass']
c_Xe = speed_of_sound(gamma_Xe,T,M_Xe)

# I
M_I = 0.12690447 # kg/mol
m_I = available_particles['I']['mass']
c_I = speed_of_sound(gamma_Xe,T,M_Xe)

# general params
drift = 30 # m/s
Sin = 0.001*0.001 # m2
n_out = 1.75e19 # m-3, output

def flux(n, Sin, T, m, drift):
    v_mean = get_maxwellian_mean_speed_from_temperature(T,m)
    return n*Sin*(0.25*v_mean+drift)

flux_out = flux(n_out, Sin, T, m_Xe, drift)


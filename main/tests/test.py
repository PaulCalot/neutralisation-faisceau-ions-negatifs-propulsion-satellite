import scipy.integrate as integrate
import numpy as np
from main import CONSTANTS
from main import MyVector

k = CONSTANTS['BOLTZMAN_CONSTANT']
    
def speed_density(drift, speed, m, temperature):
    vx = speed.x-drift.x
    vy = speed.y-drift.y
    vz = speed.z-drift.z
    factor = -m/(2*k*temperature)
    return np.exp(-factor*(vx*vx+vy*vy+vz*vz))

def normalized_speed_density(drift, speed, m, temperature):
    density = speed_density(drift, speed, m, temperature)
    factor = (m/(2*np.pi*k*temperature))**(3/2)
    return factor*density

def speed_distribution(drift, speed, m, temperature):
    pre_factor = (m/(2*np.pi*k*temperature))**(3/2)
    integrate.quad(lambda vx,vy,vz : speed_density(drift, MyVector(vx,vy,vz), m, temperature))
    
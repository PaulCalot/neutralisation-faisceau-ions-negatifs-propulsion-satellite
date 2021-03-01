import numpy as np
import scipy.stats as ss
import warnings
from math import pi as math_pi
from dolfin import Point
from main import cfg_tools
# TODO : make it better (a class ?? -> to get in an easier manner all the things we may need)
def segment(Point1,Point2):
    """
    Renvoie une liste des 4 coordonnées des extrémités à partir des deux points donnés
    [x1,y1,x2,y2]
    """
    return [Point1[0],Point1[1],Point2[0],Point2[1]]

ATOMIC_MASS = 1.66053906660e-27 # kg
NUCLEON_MASS = 1.672e-27 # kg
ELECTRON_MASS = 9.11e-31
BOLTZMAN_CONSTANT = 1.38064e-23 # J K−1
ELECTRON_CHARGE = -1.6e-19 # C
ELECTRON_EFFECTIVE_DIAMETER = 2.8179403227e-15 # m
CONSTANTS = {
    'ATOMIC_MASS':ATOMIC_MASS,
    'NUCLEON_MASS':NUCLEON_MASS,
    'ELECTRON_MASS':ELECTRON_MASS,
    'BOLTZMAN_CONSTANT':BOLTZMAN_CONSTANT,
    'ELECTRON_CHARGE':ELECTRON_CHARGE,
    'ELECTRON_EFFECTIVE_DIAMETER':ELECTRON_EFFECTIVE_DIAMETER
}

def get_mass_part(electrons_nb, protons_number, neutrons_number):
    return (neutrons_number+protons_number)*NUCLEON_MASS+electrons_nb*ELECTRON_MASS

available_particles = {
    'I':{
        'symbol' : 'I',
        'name' : 'iodine',
        'atomic mass': 126.90447,
        'mass' : get_mass_part(53, 53, 74),
        'charge' : 0,
        'effective diameter' : 4e-10
    },
    'I+':{
        'symbol' : 'I+',
        'name' : 'iodine cation',
        'mass' : get_mass_part(52, 53, 74),
        'charge' : -ELECTRON_CHARGE,
        'effective diameter' : 4e-10
    },
    'I-':{
        'symbol' : 'I-',
        'name' : 'iodide',
        'mass' : get_mass_part(54, 53, 74),
        'charge' : ELECTRON_CHARGE,
        'effective diameter' : 4e-10
    },
    'I2':{
        'symbol' : 'I2',
        'name' : 'diiodine',
        'mass' : 2*get_mass_part(53, 53, 74),
        'charge' : 0,
        # TODO
        'effective diameter' : 4e-10 + 266e-12 # distance between the two atoms
    },
    'Xe':{
        'symbol' : 'Xe',
        'name' : 'xenon',
        'atomic mass': 131.293,
        'mass' : get_mass_part(54, 54, 74),
        'charge' : 0,
        'effective diameter' : 4.32e-10 # Van der Waals radius = 216 pm
    },
    'e':{
        'symbol' : 'e',
        'name' : 'electron',
        'mass' : ELECTRON_MASS,
        'charge' : ELECTRON_CHARGE,
        'effective diameter' : ELECTRON_EFFECTIVE_DIAMETER
    }
}

# return the parameters use in scipy.stats.maxwell
def get_maxwellian_params(μ=0, σ=1):
    a = σ * np.sqrt(np.pi/(3.0*np.pi - 8.0)) # https://mathworld.wolfram.com/MaxwellDistribution.html
    m = 2.0*a*np.sqrt(2.0/np.pi)
    loc = μ - m
    return loc, a

def get_gaussian_params_maxwellian(T,m):
    #a = np.sqrt(BOLTZMAN_CONSTANT*T/m)
    #mu = 2*a*np.sqrt(2.0/np.pi)
    #sigma = a*a*(3*np.pi-8)/np.pi
    v_mean = np.sqrt(3*BOLTZMAN_CONSTANT*T/m)
    return v_mean

def get_maxwellian_mean_speed_from_temperature(T,m):
    return np.sqrt(8.0*BOLTZMAN_CONSTANT*T/(np.pi*m))

def sort_segments(walls):
    # sorting walls by x and then y
    segments = []
    for i, wall in enumerate(walls):
        x1,y1,x2,y2 = wall
        min_x, max_x, min_y, max_y = x1, x2, y1, y2
    
        if(x1 > x2):
            min_x, max_x = x2, x1
            min_y, max_y = y2, y1
        elif(x1 < x2):
            min_x, max_x = x1, x2
            min_y, max_y = y1, y2
        else :
            if(y1 > y2):
                min_x, max_x = x2, x1
                min_y, max_y = y2, y1
            elif(y2 > y1):
                min_x, max_x = x1, x2
                min_y, max_y = y1, y2
            else :
                warnings.warn("The {}th segment-wall has the same two extremities : (x1,y1) = ({},{}) and (x2,y2) = ({},{}).".format(i,x1,y1,x2,y2))
        segments.append(segment(Point(min_x, min_y),Point(max_x, max_y)))
    return segments

def collision_frequency_th_T(T, m, d, n):
    return 0.5*n*4*np.sqrt(BOLTZMAN_CONSTANT*T/m)*d*d*n

def collision_frequency_th_V(v_mean, mean_free_path, n):
    return 0.5*n*v_mean/mean_free_path

def collision_frequency_exp(dt, number_of_dt, number_of_collisions):
    return number_of_collisions/(dt*number_of_dt)


def mean_free_path(effective_diameter, particle_density):
    mean_free_path = 1/(np.sqrt(2)*np.pi*effective_diameter*effective_diameter*particle_density)
    return mean_free_path

def get_min_mean_free_path(effective_diameters, particles_densities):
    min_mean_free_path = 100000
    for effective_diameter, density in zip(effective_diameters, particles_densities):
        min_mean_free_path = min(min_mean_free_path, mean_free_path(effective_diameter, density))
    return min_mean_free_path

# def convert_to_string_():
#     if(type==type(1)):
#         return str(int)

def convert_list_to_string(mylist):
    newlist = [str(item) for item in mylist]
    return ','.join(newlist)

def convert_string_to_list(mystring, sep = ','): # by default ',' sep
    if(type(mystring) != str): 
        if(type(mystring) == float or int) : return [mystring]
        else : 
            print('Got unexpected argument : {} of type {}'.format(mystring, type(mystring)))
            raise ValueError
    return list(map(cfg_tools.read_args_multiple_types, mystring.split(sep)))

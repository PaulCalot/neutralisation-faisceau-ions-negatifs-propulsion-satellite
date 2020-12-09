import numpy as np
import scipy.stats as ss
import warnings
from math import pi as math_pi
from dolfin import Point

# TODO : make it better (a class ?? -> to get in an easier manner all the things we may need)
def segment(Point1,Point2):
    """
    Renvoie une liste des 4 coordonnées des extrémités à partir des deux points donnés
    [x1,y1,x2,y2]
    """
    return [Point1[0],Point1[1],Point2[0],Point2[1]]


NUCLEON_MASS = 1.672e-27 # kg
ELECTRON_MASS = 9.11e-31
def get_mass_part(electrons_nb, protons_number, neutrons_number):
    return (neutrons_number+protons_number)*NUCLEON_MASS+electrons_nb*ELECTRON_MASS


# return the parameters use in scipy.stats.maxwell
def get_maxwellian_params(μ=0, σ=1):
    a = σ * np.sqrt(np.pi/(3.0*np.pi - 8.0)) # https://mathworld.wolfram.com/MaxwellDistribution.html
    m = 2.0*a*np.sqrt(2.0/np.pi)
    loc = μ - m
    return loc, a

# not sorting correctly yet
def sort_segments(walls):
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

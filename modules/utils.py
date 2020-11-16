import numpy as np
import scipy.stats as ss

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

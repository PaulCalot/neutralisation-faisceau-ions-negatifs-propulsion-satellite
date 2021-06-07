from __future__ import print_function

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import mshr
import numpy as np
import scipy.integrate as integrate
from fenics import *
import pickle

# local imports
from modules.mesh_utils import get_mesh
from modules.physics_utils import get_VandE, compute_trajectory, coord_impact, intersection, compute_trajectory_1dt


mesh_dict = { # we need L < L_mot
    'L_mot' : .005,
    'l_mot' : .003,
    'L_1' : .003, # dim_trou = L_mot-L = .002 ici
    'l_1': .001,
    'L_2' : .004, # dim_trou = .001 ici
    'l_2' : .01,
    'delta_vert_12' : .001,
    'L_vacuum' : .05,
    'l_vacuum': .02,
    'mesh_resolution' : 100,
    'refine_mesh' : True,
}

phi_dict = {
    'Phi_top_mot' : 0,
    'Phi_bord_mot': 'N',
    'Phi_electrode1' :30,
    'Phi_inter_electrode':'N',
    'Phi_electrode2':300,
    'Phi_sup_vacuum':'N',
    'Phi_inf_vacuum':'N',
}

physics_consts_dict = {
    'rhoelec': 0,
    'l_rho':0, #dist between rho and 0 // <l_2 // if 0, we consider rhoelec uniform
    'PERMITTIVITY' : 8.54e-12,
    'CHARGE':1.6e-19,#Best ==> 'tmax' : .00002 and 'dt' : .00000001 and '%Nout':100
    'M_NUCLEON':1.7e-27,
}
#Best ==> 'tmax' : .00002 and 'dt' : .00000001 and '%Nout':100
integration_parameters_dict = { 
    'tmax' : .00002,
    'dt' : .00000001,
    '%Nout': 100,
}

injection_dict = {
    'Nombre de particules':10, # 10 to start with
    'proportion de I':0,
    'proportion de I+':0,
    'proportion de I-':1,
    'débit de particule en entrée de la grille':1e9,
    'gamma':1,
    'Vion':2000,
    'Vneutre':200,
    'Sigma_Vion':500,
    'Sigma_Vneutre':80,
}

mode_dict={
    'Choc symétrique?':False,
    'Contact inter particules?':False,
    'coef inelasticite':0.1,
    'coef scattering':0.3,
    'proba perte q par choc':1,
    'X0 fixe ?':False,
}

mesh, segments_list, zone = get_mesh(mesh_dict)
Phi, E, f = get_VandE(mesh, mesh_dict, phi_dict, physics_consts_dict)
Ex, Ey = E.split(deepcopy=True)
NE=sqrt(dot(E,E))

Nb_out, pf_3_especes, alpha_moy_3_especes, alpha_sigma_3_especes, V_moy_3_especes, V_sigma_3_especes, \
listes_x, listes_y, listes_vx, listes_vy, listes_q, liste_t, \
liste_alpha1, liste_alpha2, liste_alpha3, liste_V1, liste_V2, liste_V3 = \
compute_trajectory(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list, zone, E, physics_consts_dict, True)

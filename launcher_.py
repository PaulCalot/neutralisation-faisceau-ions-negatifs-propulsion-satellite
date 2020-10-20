from __future__ import print_function

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import mshr
import numpy as np
import scipy.integrate as integrate
from fenics import *

from modules.mesh_utils import get_mesh
from modules.physics_utils_pas_modif import get_VandE, compute_trajectory_bigN, compute_trajectory_smallN

# defining dicts, we need L_mot-L_1 < 2*l_mot !
mesh_dict = {
    'L_mot' : .01,
    'l_mot' : .003,
    'L_1' : .0045,
    'l_1': .005,
    'L_2' : .007,
    'l_2' : .015,
    'delta_vert_12' : .005,
    'L_vacuum' : .1,
    'l_vacuum': .05,
    'mesh_resolution' : 100,
    'refine_mesh' : True,
}

phi_dict = {
    'Phi_top_mot' : 0,
    'Phi_bord_mot': 'N',
    'Phi_electrode1' :100,
    'Phi_inter_electrode':'N',
    'Phi_electrode2':300,
    'Phi_sup_vacuum':'N',
    'Phi_inf_vacuum':'N',
}

physics_consts_dict = {
    'rhoelec': 0,
    'PERMITTIVITY' : 8.54e-12,
}

mesh, segments_list, zone = get_mesh(mesh_dict)

Phi, E = get_VandE(mesh, mesh_dict, phi_dict, physics_consts_dict)
Ex, Ey = E.split(deepcopy=True)
NE=sqrt(dot(E,E))

plt.figure(figsize=(10,10))
plot(Phi)
plt.show()

integration_parameters_dict = {
    'tmax' : .00001,
    'dt' : .0000001,
}


injection_dict = {
    'Nombre de particules':10,
    'proportion de I':0,
    'proportion de I+':0,
    'proportion de I-':1,
    'débit de particule en entrée de la grille':1e9,
}

mode_dict={
    'Elastique?':True,
    'Transfert de charge?':True,
    'Contact inter particules?':False,
    'perte u par contact':0.05,
    'proba perte q par contact':0.4,
}

from time import time

nb_simulations = 10
time_smallN = np.zeros(nb_simulations)
time_bigN = np.zeros(nb_simulations)

for k in range(nb_simulations):
    t = time()
    liste_pf, liste_alpha, liste_V, listes_x, listes_y, listes_vx, listes_vy, listes_q, liste_t = \
    compute_trajectory_smallN(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list, zone, E)
    time_smallN[k] = time()-t
    
    
for k in range(nb_simulations):
    t = time()
    liste_pf, liste_alpha, liste_V = \
    compute_trajectory_bigN(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list, zone,E)
    time_bigN[k] = time()-t
    
print("Mean time for {} simulations and {} particules : \n\t for the smallN algorithm : {} seconds \n\t\t versus \n\t for the bigN algorithm : {} seconds".format(nb_simulations,injection_dict['Nombre de particules'],round(time_smallN.mean(),1),round(time_bigN.mean(),1)))
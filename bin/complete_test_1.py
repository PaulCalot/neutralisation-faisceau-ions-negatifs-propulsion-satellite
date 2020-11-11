# TODO : nettoyer le code ; ajouter params en entrée lorsqu'on lance (ou via un cfg plus tard)

from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

from fenics import plot # to plot the mesh
from fenics import sqrt, dot

# local imports
from modules.mesh_utils import get_mesh
from modules.physics_utils import get_VandE, compute_trajectory
#from modules.plotting_utils import ?

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

integration_parameters_dict = {
    'tmax' : 4e-6,
    'dt' : 5e-10,
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
    'Contact inter particules?':True,
    'perte u par contact':0.00,
    'proba perte q par contact':0.4,
}

liste_pf, liste_alpha, liste_V, listes_x, listes_y, listes_vx, listes_vy, listes_q, liste_t = \
    compute_trajectory(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list,
    zone, E, save_trajectory = True, verbose = False,  animate = True)

plt.figure(figsize=(20,20))
fig=plot(NE)
for i in range (injection_dict['Nombre de particules']):
    plt.plot(listes_x[i],listes_y[i],linestyle='-',color='r')
    plt.scatter(listes_x[i][0],listes_y[i][0],color='r')
    plt.scatter(listes_x[i][-1],listes_y[i][-1],color='r')
plt.title('E (V/m)', size=35)
plt.xlabel('x (m)',size=10)
plt.ylabel('y (m)',size=10)
fig.set_cmap("viridis")
plt.colorbar(fig)
plt.savefig("test_1.png", dpi=400)

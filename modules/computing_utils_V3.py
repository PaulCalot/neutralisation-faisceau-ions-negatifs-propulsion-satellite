from __future__ import print_function

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import dolfin
import mshr
import numpy as np
from modules.mesh_utils_V2 import get_mesh
from modules.physics_utils_V3 import compute_trajectory,get_VandE
from fenics import *
import pickle

def main(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, phi_dict, physics_consts_dict, liste_couples, first_time=False):
    
    if first_time==True:
        liste_couples_testes = []
        liste_Nb_out=[]
        liste_pf = []
        liste_alpha_moy = []
        liste_alpha_sigma = []
        liste_V_moy = []
        liste_V_sigma = []
        liste_alpha_detail = []
        liste_V_detail = []
        liste_V_detail_inside = []
        
    else:
        with open("CAS_TEST_V2_liste_couples_testes.txt", "rb") as fp:
            liste_couples_testes = pickle.load(fp)
        with open("CAS_TEST_V2_liste_Nb_out.txt", "rb") as fp:
            liste_Nb_out = pickle.load(fp)
        with open("CAS_TEST_V2_liste_pf.txt", "rb") as fp:
            liste_pf = pickle.load(fp)
        with open("CAS_TEST_V2_liste_alpha_moy.txt", "rb") as fp:
            liste_alpha_moy = pickle.load(fp)
        with open("CAS_TEST_V2_liste_alpha_sigma.txt", "rb") as fp:
            liste_alpha_sigma = pickle.load(fp)
        with open("CAS_TEST_V2_liste_V_moy.txt", "rb") as fp:
            liste_V_moy = pickle.load(fp)
        with open("CAS_TEST_V2_liste_V_sigma.txt", "rb") as fp:
            liste_V_sigma = pickle.load(fp)
        with open("CAS_TEST_V2_liste_alpha_detail.txt", "rb") as fp:
            liste_alpha_detail = pickle.load(fp)
        with open("CAS_TEST_V2_liste_V_detail.txt", "rb") as fp:
            liste_V_detail = pickle.load(fp)
        with open("CAS_TEST_V2_liste_V_detail_inside.txt", "rb") as fp:
            liste_V_detail_inside = pickle.load(fp)

##############################################################################################################################     
    for couple in liste_couples:

        if couple in liste_couples_testes:
            print('combinaison '+str(couple)+' déjà testée')

        else:    
            mesh_dict['L_2'] = couple[0]
            mesh_dict['l_2'] = couple[1]

            mesh, segments_list, zone = get_mesh(mesh_dict)
            Phi, E, f = get_VandE(mesh, mesh_dict, phi_dict, physics_consts_dict)

            Nb_out, pf_2_especes, alpha_moy_2_especes, alpha_sigma_2_especes, V_moy_2_especes, V_sigma_2_especes, \
listes_x, listes_y, listes_vx, listes_vy, listes_q, liste_t, \
liste_alpha_neutre, liste_alpha_ion, liste_V_neutre, liste_V_ion,\
liste_V_neutre_inside, liste_V_ion_inside= \
compute_trajectory(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list, zone, E, physics_consts_dict, True, True)

            liste_Nb_out.append(Nb_out)
            liste_pf.append(pf_2_especes)
            liste_alpha_moy.append(alpha_moy_2_especes)
            liste_alpha_sigma.append(alpha_sigma_2_especes)
            liste_V_moy.append(V_moy_2_especes)
            liste_V_sigma.append(V_sigma_2_especes)
            liste_couples_testes.append(couple)
            liste_alpha_detail.append([liste_alpha_neutre, liste_alpha_ion])
            liste_V_detail.append([liste_V_neutre, liste_V_ion])
            liste_V_detail_inside.append([liste_V_neutre_inside, liste_V_ion_inside])


            print('combinaison '+str(couple)+' testée')
            print('------------------------------------------------------')

            with open("CAS_TEST_V2_liste_couples_testes.txt", "wb") as fp:
                pickle.dump(liste_couples_testes, fp)
            with open("CAS_TEST_V2_liste_Nb_out.txt", "wb") as fp:
                pickle.dump(liste_Nb_out, fp)
            with open("CAS_TEST_V2_liste_pf.txt", "wb") as fp:
                pickle.dump(liste_pf, fp)
            with open("CAS_TEST_V2_liste_alpha_moy.txt", "wb") as fp:
                pickle.dump(liste_alpha_moy, fp)
            with open("CAS_TEST_V2_liste_alpha_sigma.txt", "wb") as fp:
                pickle.dump(liste_alpha_sigma, fp)
            with open("CAS_TEST_V2_liste_V_moy.txt", "wb") as fp:
                pickle.dump(liste_V_moy, fp)
            with open("CAS_TEST_V2_liste_V_sigma.txt", "wb") as fp:
                pickle.dump(liste_V_sigma, fp)
            with open("CAS_TEST_V2_liste_alpha_detail.txt", "wb") as fp:
                pickle.dump(liste_alpha_detail, fp)
            with open("CAS_TEST_V2_liste_V_detail.txt", "wb") as fp:
                pickle.dump(liste_V_detail, fp)
            with open("CAS_TEST_V2_liste_V_detail_inside.txt", "wb") as fp:
                pickle.dump(liste_V_detail_inside, fp)

    
    
    
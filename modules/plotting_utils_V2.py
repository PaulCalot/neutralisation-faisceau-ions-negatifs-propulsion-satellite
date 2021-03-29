from __future__ import print_function

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import dolfin
import mshr
import numpy as np
from modules.mesh_utils_V2 import get_mesh
from fenics import *
# TODO : delete unused imporrt

def plot_elements(mesh,E,Phi,mesh_dict,quality):
    
    Ex, Ey = E.split(deepcopy=True)
    NE=sqrt(dot(E,E))
    
    borne_inf_x=-1.5*mesh_dict['L_mot']
    borne_sup_x=1.5*mesh_dict['L_mot']
    borne_inf_y=-2.5*(.5*mesh_dict['l_mot']+mesh_dict['l_1']+mesh_dict['l_2']\
                      +mesh_dict['delta_vert_12']+mesh_dict['l_3']+mesh_dict['delta_vert_23'])
    borne_sup_y=.5*mesh_dict['l_mot']

    if quality==False:
        
        plt.figure(figsize=(20,10))

        plt.subplot(1,3,1)
        fig=plot(mesh)
        plt.title('mesh', size=30)
        plt.xlim(borne_inf_x,borne_sup_x)
        plt.ylim(borne_inf_y,borne_sup_y)
        plt.xlabel('x(m)', size=30)
        plt.ylabel('y(m)', size=30)

        plt.subplot(1,3,2)
        fig=plot(Phi)
        plt.title('U(V)', size=30)
        fig.set_cmap("viridis") 
        plt.colorbar(fig)
        plt.xlim(borne_inf_x,borne_sup_x)
        plt.ylim(borne_inf_y,borne_sup_y)
        plt.xlabel('x(m)', size=30)
        plt.ylabel('y(m)', size=30)

        plt.subplot(1,3,3)
        fig=plot(NE)
        plt.title('E(V/m)', size=30)
        fig.set_cmap("viridis") 
        plt.colorbar(fig)
        plt.xlim(borne_inf_x,borne_sup_x)
        plt.ylim(borne_inf_y,borne_sup_y)
        plt.xlabel('x(m)', size=30)
        plt.ylabel('y(m)', size=30)

        plt.show()
    
    if quality==True:
        
        plt.figure(figsize=(20,10))
        fig=plot(mesh)
        plt.title('mesh', size=30)
        plt.xlabel('x(m)', size=30)
        plt.ylabel('y(m)', size=30)
        plt.show()
        
        plt.figure(figsize=(20,10))

        plt.subplot(1,2,1)
        fig=plot(Phi)
        plt.title('U(V)', size=30)
        fig.set_cmap("viridis") 
        plt.colorbar(fig)
        plt.xlim(borne_inf_x,borne_sup_x)
        plt.ylim(borne_inf_y,borne_sup_y)
        plt.xlabel('x(m)', size=30)
        plt.ylabel('y(m)', size=30)

        plt.subplot(1,2,2)
        fig=plot(NE)
        plt.title('E(V/m)', size=30)
        fig.set_cmap("viridis") 
        plt.colorbar(fig)
        plt.xlim(borne_inf_x,borne_sup_x)
        plt.ylim(borne_inf_y,borne_sup_y)
        plt.xlabel('x(m)', size=30)
        plt.ylabel('y(m)', size=30)

        plt.show()
        

def basic_plot_trajectories(E,listes_x,listes_y):
    NE=sqrt(dot(E,E))
    plt.figure(figsize=(12,12))
    fig=plot(NE)
    for i in range (len(listes_x)):
        plt.plot(listes_x[i],listes_y[i],linestyle='-',color='r')
        plt.scatter(listes_x[i][0],listes_y[i][0],color='r')
        plt.scatter(listes_x[i][-1],listes_y[i][-1],color='r')
    plt.title('E (V/m)', size=35)
    plt.xlabel('x (m)',size=35)
    plt.ylabel('y (m)',size=35)
    fig.set_cmap("viridis")
    plt.colorbar(fig)
    plt.show()
        

def plot_cas_test(mesh_dict, liste_couples):
    for couple in liste_couples: 
        mesh_dict['L_2']=couple[0]
        mesh_dict['l_2']=couple[1]
        mesh, segments_list, zone = get_mesh(mesh_dict)
        plt.figure(figsize=(10,10))
        if couple[1]==0.001:
            plt.ylim(-0.025,0.003)
        else:
            plt.ylim(-0.035,0.003)
        plot(mesh)
        plt.xlabel('x(m)', size=30)
        plt.ylabel('y(m)', size=30)
        plt.show()
        
def neutralisation_cas_test(liste_pf):
    abscisse=[1,2,3,4]
    ordonnee=liste_pf[:,0]

    plt.figure(figsize=(10,10))
    plt.scatter(abscisse, ordonnee, color='k')
    plt.xlabel('N° du cas test', size=30)
    plt.ylabel("taux de neutralisation en sortie",size=30)
    plt.ylim(-0.01,1)
    plt.xticks([1, 2, 3, 4])

    plt.show()

def vitesse_par_espece_cas_test(liste_V_moy,liste_V_sigma,liste_pf,liste_Nb_out):
    abscisse=[1,2,3,4]
    ordonnee1=liste_V_moy[:,0]
    barre1=liste_V_sigma[:,0]
    barre1=1.96*barre1*np.power(liste_pf[:,0]*liste_Nb_out,-.5)
    ordonnee2=liste_V_moy[:,1]
    barre2=liste_V_sigma[:,1]
    barre2=1.96*barre2*np.power(liste_pf[:,1]*liste_Nb_out,-.5)

    plt.figure(figsize=(20,7))

    plt.subplot(1,2,1)
    plt.scatter(abscisse, ordonnee1, color='k')
    plt.errorbar(abscisse, ordonnee1, yerr=barre1, fmt = 'none', capsize = 7, elinewidth=1, ecolor='k')
    plt.xlabel('N° du cas test', size=30)
    plt.ylabel('vitesse moyenne (m/s)',size=30)
    plt.title('I',size=30)
    plt.ylim(13000,20000)
    plt.xticks([1, 2, 3, 4])

    plt.subplot(1,2,2)
    plt.scatter(abscisse, ordonnee2, color='k')
    plt.errorbar(abscisse, ordonnee2, yerr=barre2, fmt = 'none', capsize = 7, elinewidth=1, ecolor='k')
    plt.xlabel('N° du cas test', size=30)
    plt.ylabel('vitesse moyenne (m/s)',size=30)
    plt.title('I-',size=30)
    plt.ylim(13000,20000)
    plt.xticks([1, 2, 3, 4])
    
    plt.show()

def angle_par_espece_cas_test(liste_alpha_moy,liste_alpha_sigma,liste_pf,liste_Nb_out):
    abscisse=[1,2,3,4]
    ordonnee1=liste_alpha_moy[:,0]
    barre1=liste_alpha_sigma[:,0]
    barre1=1.96*barre1*np.power(liste_pf[:,0]*liste_Nb_out,-.5)
    ordonnee2=liste_alpha_moy[:,1]
    barre2=liste_alpha_sigma[:,1]
    barre2=1.96*barre2*np.power(liste_pf[:,1]*liste_Nb_out,-.5)

    plt.figure(figsize=(20,7))

    ax = plt.subplot(1,2,1)
    plt.scatter(abscisse, ordonnee1, color='k')
    plt.errorbar(abscisse, ordonnee1, yerr=barre1, fmt = 'none', capsize = 7, elinewidth=1, ecolor='k')
    y_pi   = ordonnee1/np.pi
    unit   = 0.05
    y_tick = np.arange(-0.05, 0.051, unit)
    y_label = [r"$-\frac{\pi}{20}$", r"$0$", r"$+\frac{\pi}{20}$"]
    ax.set_yticks(y_tick*np.pi)
    ax.set_yticklabels(y_label, fontsize=20)
    plt.xlabel('N° du cas test', size=30)
    plt.ylabel('angle (rad) du flux',size=30)
    plt.ylim(-np.pi/15, np.pi/15)
    plt.xticks([1, 2, 3, 4])
    plt.title('I',size=30)

    ax = plt.subplot(1,2,2)
    plt.scatter(abscisse, ordonnee2, color='k')
    plt.errorbar(abscisse, ordonnee2, yerr=barre2, fmt = 'none', capsize = 7, elinewidth=1, ecolor='k')
    y_pi   = ordonnee2/np.pi
    unit   = 0.05
    y_tick = np.arange(-0.05, 0.051, unit)
    y_label = [r"$-\frac{\pi}{20}$", r"$0$", r"$+\frac{\pi}{20}$"]
    ax.set_yticks(y_tick*np.pi)
    ax.set_yticklabels(y_label, fontsize=20)
    plt.xlabel('N° du cas test', size=30)
    plt.ylabel("angle (rad) du flux",size=30)
    plt.ylim(-np.pi/15, np.pi/15)
    plt.xticks([1, 2, 3, 4])
    plt.title('I-',size=30)

    plt.show()
    
def vitesse_total_cas_test(liste_V_moy,liste_V_sigma,liste_pf,liste_Nb_out):
    abscisse=[1,2,3,4]
    ordonnee=liste_V_moy[:,0]*liste_pf[:,0]+liste_V_moy[:,1]*liste_pf[:,1]
    barre=liste_V_sigma[:,0]*liste_pf[:,0]+liste_V_sigma[:,1]*liste_pf[:,1]
    barre=1.96*barre*np.power(liste_Nb_out,-.5)

    plt.figure(figsize=(14,14))
    plt.suptitle('Flux Total',size=30)
    plt.scatter(abscisse, ordonnee, color='k')
    plt.errorbar(abscisse, ordonnee, yerr=barre, fmt = 'none', capsize = 7, elinewidth=1, ecolor='k')
    plt.xlabel('N° du cas test', size=30)
    plt.ylabel('vitesse moyenne (m/s)',size=30)
    plt.xticks([1, 2, 3, 4])
    plt.show()

def angle_total_cas_test(liste_alpha_moy,liste_alpha_sigma,liste_pf,liste_Nb_out):
    abscisse=[1,2,3,4]
    ordonnee=liste_alpha_moy[:,0]*liste_pf[:,0]+liste_alpha_moy[:,1]*liste_pf[:,1]
    barre=liste_alpha_sigma[:,0]*liste_pf[:,0]+liste_alpha_sigma[:,1]*liste_pf[:,1]
    barre=1.96*barre*np.power(liste_Nb_out,-.5)
    
    plt.figure(figsize=(14,14))
    ax = plt.subplot(1,1,1)
    plt.suptitle('Flux Total',size=30)
    plt.scatter(abscisse, ordonnee, color='k')
    plt.errorbar(abscisse, ordonnee, yerr=barre, fmt = 'none', capsize = 7, elinewidth=1, ecolor='k')
    y_pi   = ordonnee/np.pi
    unit   = 0.05
    y_tick = np.arange(-0.05, 0.051, unit)
    y_label = [r"$-\frac{\pi}{20}$", r"$0$", r"$+\frac{\pi}{20}$"]
    ax.set_yticks(y_tick*np.pi)
    ax.set_yticklabels(y_label, fontsize=20)
    plt.xlabel('N° du cas test', size=30)
    plt.ylabel('angle (rad) du flux',size=30)
    plt.xticks([1, 2, 3, 4])
    plt.show()

def distribution_vitesse_1_cas(liste_V_detail,n):
    plt.figure(figsize=(20,7))
    
    plt.suptitle("Distributions de vitesse pour le cas test "+str(n),size=30)
    
    plt.subplot(1,2,1)
    ordonnee=liste_V_detail[n-1,0]
    plt.hist(ordonnee,bins=15)
    plt.xlabel('vitesse moyenne (m/s)',size=20)
    plt.ylabel('effectif',size=20)
    plt.xlim(3000,23000)
    plt.title('I',size=20)
    
    plt.subplot(1,2,2)
    ordonnee=liste_V_detail[n-1,1]
    plt.hist(ordonnee,bins=15)
    plt.xlabel('vitesse moyenne (m/s)',size=20)
    plt.ylabel('effectif',size=20)
    plt.xlim(3000,23000)
    plt.title('I-',size=20)
    
    plt.show()

    
    
def distribution_angle_1_cas(liste_alpha_detail,n):
    plt.figure(figsize=(20,7))
    
    plt.suptitle("Distributions d'angle pour le cas test "+str(n),size=30)
    
    ax = plt.subplot(1,2,1)
    ordonnee=liste_alpha_detail[n-1,0]
    plt.hist(ordonnee)
    x_pi   = ordonnee/np.pi
    unit   = 0.1
    x_tick = np.arange(-0.4, 0.401, unit)
    x_label = [r"$-\frac{2\pi}{5}$", r"$-\frac{3\pi}{10}$",r"$-\frac{\pi}{5}$",r"$-\frac{\pi}{10}$", r"$0$", r"$+\frac{\pi}{10}$",r"$+\frac{\pi}{5}$",r"$+\frac{3\pi}{10}$",r"$+\frac{2\pi}{5}$"]
    ax.set_xticks(x_tick*np.pi)
    ax.set_xticklabels(x_label, fontsize=20)
    plt.xlabel('angle (rad) du flux',size=20)
    plt.ylabel('effectif',size=20)
    plt.title('I',size=20)
    
    ax = plt.subplot(1,2,2)
    ordonnee=liste_alpha_detail[n-1,1]
    plt.hist(ordonnee)
    x_pi   = ordonnee/np.pi
    unit   = 0.1
    x_tick = np.arange(-0.4, 0.401, unit)
    x_label = [r"$-\frac{2\pi}{5}$", r"$-\frac{3\pi}{10}$",r"$-\frac{\pi}{5}$",r"$-\frac{\pi}{10}$", r"$0$", r"$+\frac{\pi}{10}$",r"$+\frac{\pi}{5}$",r"$+\frac{3\pi}{10}$",r"$+\frac{2\pi}{5}$"]
    ax.set_xticks(x_tick*np.pi)
    ax.set_xticklabels(x_label, fontsize=20)
    plt.xlabel('angle (rad) du flux',size=20)
    plt.ylabel('effectif',size=20)
    plt.title('I-',size=20)
    
    plt.show()
    
def distribution_vitesse_1_cas_inside(liste_V_detail_inside,n):
    plt.figure(figsize=(20,7))
    
    plt.suptitle("Distributions de vitesse pour le cas test "+str(n),size=30)
    
    plt.subplot(1,2,1)
    ordonnee=liste_V_detail_inside[n-1,0]
    plt.hist(ordonnee,bins=15,color='orange')
    plt.xlabel('vitesse moyenne (m/s)',size=20)
    plt.ylabel('effectif',size=20)
    plt.xlim(3000,23000)
    plt.title('I',size=20)
    
    plt.subplot(1,2,2)
    ordonnee=liste_V_detail_inside[n-1,1]
    plt.hist(ordonnee,bins=15,color='orange')
    plt.xlabel('vitesse moyenne (m/s)',size=20)
    plt.ylabel('effectif',size=20)
    plt.xlim(3000,23000)
    plt.title('I-',size=20)
    
    plt.show()
    
from __future__ import print_function

# local import

from .particules import Particule
from .vector import MyVector
from .collisions_handler import CollisionHandler

# animation
import matplotlib.pyplot as plt
import matplotlib.animation as animation


import mshr
import numpy as np
from fenics import *

import sys # to get the max float value
from time import time # to time out if need
from tqdm import tqdm
u=1.7e-27
e=1.6e-19

# ------------------------- E and V computation ---------------------------- #
def get_VandE(mesh, mesh_dict, phi_dict, physics_consts_dict):
    
    V = FunctionSpace(mesh, 'P', 1)
    W = VectorFunctionSpace(mesh, 'P', 1)
    
    L_mot = mesh_dict['L_mot']
    l_mot = mesh_dict['l_mot']
    L_vacuum = mesh_dict['L_vacuum']
    l_vacuum = mesh_dict['l_vacuum']
    L_1 = mesh_dict['L_1']
    l_1 = mesh_dict['l_1']
    L_2 = mesh_dict['L_2']
    l_2 = mesh_dict['l_2']
    delta_vert_12 = mesh_dict['delta_vert_12']
    mesh_resolution = mesh_dict['mesh_resolution']
    refine_mesh = mesh_dict['refine_mesh']

    h_grid = l_1 + l_2 + delta_vert_12

    class Boundary_top_mot(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] == l_mot/2
        
    class Boundary_bord_mot(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] < l_mot/2  and x[1] > - l_mot/2
        
    class Boundary_electrode1(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] <= - l_mot/2 and x[1] >= - l_mot/2 - l_1
        
    class Boundary_inter_electrode(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] < - l_mot/2 - l_1 and x[1] > - l_mot/2 - l_1 - delta_vert_12
        
    class Boundary_electrode2(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] <= - l_mot/2 - l_1 - delta_vert_12 and x[1] >= - l_mot/2 - \
                l_1 - delta_vert_12 - l_2 and abs(x[0])<=L_mot/2
        
    class Boundary_sup_vacuum(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and (( x[1]< -l_mot/2 - h_grid and x[1] >= -l_mot/2 - \
                h_grid - l_vacuum/2) or (x[1]== -l_mot/2 - h_grid and abs(x[0])>L_mot/2))
        
    class Boundary_inf_vacuum(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1]< -l_mot/2 - h_grid - l_vacuum/2

    top_mot = Boundary_top_mot()
    bord_mot = Boundary_bord_mot()
    electrode1 = Boundary_electrode1()
    inter_electrode = Boundary_inter_electrode()
    electrode2 = Boundary_electrode2()
    sup_vacuum = Boundary_sup_vacuum()
    inf_vacuum = Boundary_inf_vacuum()

    list_Phi = [phi_dict['Phi_top_mot'],phi_dict['Phi_bord_mot'],phi_dict['Phi_electrode1'], \
        phi_dict['Phi_inter_electrode'], phi_dict['Phi_electrode2'],phi_dict['Phi_sup_vacuum'],phi_dict['Phi_inf_vacuum']]
    list_edges=[top_mot, bord_mot, electrode1, inter_electrode, electrode2, sup_vacuum, inf_vacuum]
    bc=[]
    
    for i in range(len(list_edges)):
        if list_Phi[i]!='N':
            bc.append(DirichletBC(V, Constant(list_Phi[i]), list_edges[i]))


    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(physics_consts_dict['rhoelec']/physics_consts_dict['PERMITTIVITY']) 
    a = dot(grad(u), grad(v))*dx
    L = f*v*dx
    Phi = Function(V)

    solve(a == L, Phi, bc)

    E = project(-grad(Phi), W)

    return Phi,E


# ------------------------------ Trajectory computation auxiliary functions -------------------- #

def moyenne_amelioree(liste):
    if len(liste)==0:
        return None
    return np.mean(liste)

def distrib_init(espece, mesh_dict):
    """
    Renvoie x,y,z aléatoirement en suivant une distribution propre à l'espèce
    """
    L_mot = mesh_dict['L_mot']
    l_mot = mesh_dict['l_mot']
    L_1 = mesh_dict['L_1']
    l_1 = mesh_dict['l_1']
    
    if espece=='I':
        alpha=np.random.uniform(0,2*np.pi)
        v=np.random.normal(200,1e3)
        x=np.random.uniform(-(L_mot-L_1)/2,(L_mot-L_1)/2)
        y=np.random.uniform(-l_mot/4,l_mot/4)
        return x, y, 0, v*np.cos(alpha), v*np.sin(alpha), 0
    
    elif espece=='I+':
        alpha=np.random.uniform(0,np.pi)
        v=np.random.normal(5e3,1e3)
        x=(L_mot-L_1)/2*np.cos(alpha)
        y=(L_mot-L_1)/2*np.sin(alpha)-l_mot/2
        vx=v*np.cos(alpha)
        vy=v*np.sin(alpha)
        return x, y, 0, vx, vy, 0
    
    elif espece=='I-':
        alpha=np.random.uniform(0,np.pi)
        v=np.random.normal(5e3,1e3)
        x=(L_mot-L_1)/2*np.cos(alpha)
        y=(L_mot-L_1)/2*np.sin(alpha)-l_mot/2
        vx=-v*np.cos(alpha)
        vy=-v*np.sin(alpha)
        return x, y, 0, vx, vy, 0

def f(Y,t,m,q,zone,E):
    '''
    Renvoie la dérivée de Y en y (en faisant un bilan des forces notamment) pr être entré dans RK4
    Y=[x, y, z, vx, vy, vz] ce n'est pas le Y de liste_Y
    '''
    Ex, Ey = E.split(deepcopy=True)
    vx=Y[3]
    vy=Y[4]
    vz=Y[5]
    if zone.inside(Point(Y[0],Y[1])):
        ax = (1/m) * q * Ex(Y[0], Y[1])
        ay = (1/m) * q * Ey(Y[0], Y[1])
    else :
        ax = 0  #utile si les ki st hors du mesh,
        ay = 0
    az=0
    return [vx, vy, vz, ax, ay, az]

# ------------------------------ Trajectory computation main function -------------------- #

def compute_trajectory(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list,
                       zone, E, save_trajectory = False, verbose = True, animate = True):
    """
    Renvoie la proportion d'espèce , 
    l'angle moyen, 
    et la norme de la vitesse moyenne en sortie.
    
    Ainsi que la liste des instants t et les listes de liste de x,y,vx,vy,q de chaque particule
    """
    
    tmax=integration_parameters_dict['tmax']
    dt=integration_parameters_dict['dt']
    
    l_mot=mesh_dict['l_mot']
    l_1=mesh_dict['l_1']
    l_2=mesh_dict['l_2']
    delta_vert_12=mesh_dict['delta_vert_12']
    l_vacuum=mesh_dict['l_vacuum']
    h_grid=l_1 + l_2 + delta_vert_12
    
    N=injection_dict['Nombre de particules']
    p1=injection_dict['proportion de I']
    p2=injection_dict['proportion de I+']
    p3=injection_dict['proportion de I-']
    DN=injection_dict['débit de particule en entrée de la grille']
    
    ### On initialise une liste_Y comportant N lignes, chaque ligne correspond à une particule et est un couple (Particule,Entier)
    ### l'entier vaut -1 si la particule n'est pas encore injectée, 0 si est dedans et 1 si est sortie
    
    if verbose : 
        print("INITIALIZING SIMULATION")
        print('\t max simulation time : {} seconds \n\t time step : {} seconds \n\t Number of particules : {}'.format(tmax, dt, N)) 
        print('\t particules proportions : \n\t\t I  : {:.0%} \n\t\t I+ : {:.0%} \n\t\t I- : {:.0%}'.format(p1, p2, p3))
    
    N2=int(p2*N)
    N3=int(p3*N)
    N1=N-N2-N3
    
    list_parts=[]
    for n in range (N1):
        x,y,z,vx,vy,vz=distrib_init('I', mesh_dict)
        list_parts.append(Particule(charge = 0, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I", status = 0)) # note that there still is a radius that can be set. By default it is set to IODINE_RADIUS.
    for n in range (N1,N1+N2):
        x,y,z,vx,vy,vz=distrib_init('I+', mesh_dict)
        list_parts.append(Particule(charge = e, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I+", status = 0))
    for n in range (N1+N2,N):
        x,y,z,vx,vy,vz=distrib_init('I-', mesh_dict)
        list_parts.append(Particule(charge = -e, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I-", status = 0))
    
    max_steps = int(tmax/dt)+1
    nombre_max_injecte_par_tour=int(DN*dt)+1
    #nombre_tour_plein_debit=int(N/nombre_max_injecte_par_tour)
    #nombre_derniere_injection=N-nombre_tour_plein_debit*nombre_max_injecte_par_tour

    np.random.shuffle(list_parts)
    for i in range(nombre_max_injecte_par_tour, N):
        # setting particules not yet injected to status -1
        list_parts[i].set_status(-1)
    # /!\ creating the collision handler /!\
        # we had to wait for the shuffle of liste_Y as the index we'll have is the on in liste_Y
    Cond1=mode_dict['Elastique?']
    Cond2=mode_dict['Transfert de charge?']
    Cond3=mode_dict['Contact inter particules?']
    eta=mode_dict['perte u par contact'] if Cond1 else 0 # a priori idem pr I, I+ et I-
    p=mode_dict['proba perte q par contact'] if Cond2 else 0 # a priori idem pr I+ et I-

    if(Cond3):
        collision_handler = CollisionHandler(list_parts, segments_list, f, eta, p)
    else :
        #collision_handler = CollisionHandlerWithoutInterPartCollision(list_parts, segments_list, zone, f, eta, p)
        pass
    # END
        # is it useful ?
    if(save_trajectory):
        liste_t=[0]
        listes_x=[]
        listes_y=[]
        listes_vx=[]
        listes_vy=[]
        listes_q=[]
        for n in range(N):
            particule = list_parts[n]
            pos = particule.get_pos()
            speed = particule.get_speed()
            listes_x.append([pos.x])
            listes_y.append([pos.y])
            listes_vx.append([speed.x])
            listes_vy.append([speed.y])
            listes_q.append([particule.get_charge()])
    
    if verbose : print("PARTICULES INJECTION AND SIMULATION ...", end = " ")
        
    ### On calcule les trajectoires en 2 temps
    ### On fait d'abord les tours où il y a injection au débit max en ne traitant qu'une partie des particules
    ### Puis on fait tourner en continu jusqu'à tmax ou à ce qu'il n'y ait plus de particules à l'intérieur
        
    t=0
    Nb_out=0
    
    injection_finished = False
    """ try for animating the particules
    if(animate):
        NE=sqrt(dot(E,E))
        fig = plt.figure(figsize=(10,10))
        #fig=plot(NE)
        #fig.set_cmap("viridis")
        #lt.colorbar(fig)
        ax = fig.add_subplot(111, aspect='equal')
        X_pos = []
        Y_pos = []
            
        for n in range(N):
            part = list_parts[n]
            pos = part.get_pos()

            # setting position for animating
            X_pos.append(pos.x)
            Y_pos.append(pos.y)
            
        line, = ax.plot(X_pos, Y_pos, 'ro', ms=3)
        ax.set_xlim(-0.05,0.05)
        ax.set_ylim(-0.03,0)
        
        def init():
            return line,
        
        def animate_func(i):
            print("DEBUG")
            if(Nb_out >= N):
                #return
                pass
            if (verbose): #and np.random.random_sample()<max(dt/tmax,0.05)):
                print ('elapsed time : {:.0%} , particules remaining : {:.0%}.'.format(t/tmax,1-Nb_out/N))
            if(not injection_finished):
                for n in range (i*nombre_max_injecte_par_tour, (i+1)*nombre_max_injecte_par_tour):
                    list_parts[n].set_status(0)
                    if n==N-1: 
                        if(verbose): print("INJECTION FINISHED")
                        injection_finished = True
                        #break
            collision_handler.step(dt, t, E, zone)
            collision_handler.update_all_events()
            if(verbose): print("Next time : {} sec. Time of next collision : {} sec.".format(t+dt,t+collision_handler.get_next_collision()))
            print("DEBUG"

            X_pos = []
            Y_pos = []
            
            for n in range(N):
                print("DEBUG")
                part = list_parts[n]
                pos = part.get_pos()
                
                # setting position for animating
                X_pos.append(pos.x)
                Y_pos.append(pos.y)
                
                if(save_trajectory):
                    speed = particule.get_speed()
                    listes_x[n].append(pos.x)
                    listes_y[n].append(pos.y)
                    listes_vx[n].append(speed.x)
                    listes_vy[n].append(speed.y)
                    listes_q[n].append(particule.get_charge())
                    liste_t.append(t+dt)

                if pos.y<-l_mot/2-h_grid-l_vacuum/2:
                    part.set_status(1)
                    Nb_out+=1
            t+=dt
            line.set_data(X_pos, Y_pos)
            return line,
    
        anim = animation.FuncAnimation(fig, animate_func, init_func=init, interval=dt*1000, blit=False)
        plt.show()
    """


      #for i in range(max_steps):
    for i in tqdm(range(max_steps)):
        if(Nb_out >= N):
            break
        if (verbose): #and np.random.random_sample()<max(dt/tmax,0.05)):
            print ('elapsed time : {:.0%} , particules remaining : {:.0%}.'.format(t/tmax,1-Nb_out/N))
        if(not injection_finished):
            for n in range (i*nombre_max_injecte_par_tour, (i+1)*nombre_max_injecte_par_tour):
                list_parts[n].set_status(0)
                if n==N-1: 
                    if(verbose): print("INJECTION FINISHED")
                    injection_finished = True
                    break

        collision_handler.step(dt, t, f_args = [zone, E])
        collision_handler.update_all_events()
        if(verbose): print("Next time : {} sec. Time of next collision : {} sec.".format(t+dt,t+collision_handler.get_next_collision()))
        if(injection_finished):
            for n in range(N):
                part = list_parts[n]
                pos = part.get_pos()

                if(save_trajectory):
                    speed = particule.get_speed()
                    listes_x[n].append(pos.x)
                    listes_y[n].append(pos.y)
                    listes_vx[n].append(speed.x)
                    listes_vy[n].append(speed.y)
                    listes_q[n].append(particule.get_charge())
                    liste_t.append(t+dt)

                if pos.y<-l_mot/2-h_grid-l_vacuum/2:
                    part.set_status(1)
                    Nb_out+=1
        t+=dt

    if(verbose):print('\t[OK]')

    if(verbose): 
        print("END : processing \nSTOPPING CRITERION : ", end = " ")
        if(t>tmax):
            print("TMAX REACHED ~> t = {} > {} = tmax seconds".format(t,tmax))
        elif(Nb_out==N):
            print("ALL PARTICULES EXITED ~> Nb_out = N = {}".format(Nb_out))

    ### Traitement des données de sortie
    if(verbose): print('START : data processing ...', end = " ")
    
    N1f,N2f,N3f=0,0,0
    liste_vxf1,liste_vyf1=[],[]
    liste_vxf2,liste_vyf2=[],[]
    liste_vxf3,liste_vyf3=[],[]

    for n in range(N):
        if list_parts[n].get_status()==1 and list_parts[n].get_speed().y!=0:
            particule=list_parts[n]
            speed = particule.get_speed()
            part_type = particule.get_part_type()
            if part_type=='I':
                N1f+=1
                liste_vxf1.append(speed.x)
                liste_vyf1.append(speed.y)
            elif part_type=='I+':
                N2f+=1
                liste_vxf2.append(speed.x)
                liste_vyf2.append(speed.y)
            elif part_type=='I-':
                N3f+=1
                liste_vxf3.append(speed.x)
                liste_vyf3.append(speed.y)
    
    if Nb_out!=0:
        p1f=N1f/Nb_out
        p2f=N2f/Nb_out
        p3f=N3f/Nb_out
    else:
        p1f=0
        p2f=0
        p3f=0

    # converting to arrays
    liste_vxf1=np.array(liste_vxf1)
    liste_vyf1=np.array(liste_vyf1)
    liste_vxf2=np.array(liste_vxf2)
    liste_vyf2=np.array(liste_vyf2)
    liste_vxf3=np.array(liste_vxf3)
    liste_vyf3=np.array(liste_vyf3)

    # computing quantities of interest
    liste_Vnorm1=np.power(np.power(liste_vxf1,2)+np.power(liste_vyf1,2), 1/2)
    liste_alpha1=-np.arctan(liste_vxf1/liste_vyf1)
    liste_Vnorm2=np.power(np.power(liste_vxf2,2)+np.power(liste_vyf2,2), 1/2)
    liste_alpha2=-np.arctan(liste_vxf2/liste_vyf2)
    # TODO problem here, see what Edouard made 
    liste_Vnorm3=np.power(np.power(liste_vxf3,2)+np.power(liste_vyf3,2), 1/2)
    liste_alpha3=-np.arctan(liste_vxf3/liste_vyf3)
     # ajouter z
    if(verbose): 
        print('\t[OK]')
        print('RESULTS :')
        print('\t Particules that exited : {} out of {}.'.format(Nb_out,N))
        print('\t Proportions : \n\t\t I  : {:.0%} \n\t\t I+ : {:.0%} \n\t\t I- : {:.0%}'.format(p1f,p2f,p3f))
        if(save_trajectory):
            print("=> Intermediary positions returned")
        print('\n')
        
    if(save_trajectory):
        return [p1f,p2f,p3f], \
                [moyenne_amelioree(liste_alpha1),moyenne_amelioree(liste_alpha2),moyenne_amelioree(liste_alpha3)], \
                [moyenne_amelioree(liste_Vnorm1),moyenne_amelioree(liste_Vnorm2),moyenne_amelioree(liste_Vnorm3)], \
                listes_x, listes_y, listes_vx, listes_vy, listes_q, liste_t
    return [p1f,p2f,p3f], \
            [moyenne_amelioree(liste_alpha1),moyenne_amelioree(liste_alpha2),moyenne_amelioree(liste_alpha3)], \
            [moyenne_amelioree(liste_Vnorm1),moyenne_amelioree(liste_Vnorm2),moyenne_amelioree(liste_Vnorm3)]


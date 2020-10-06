from __future__ import print_function


import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import dolfin
import mshr
import numpy as np
import scipy.integrate as integrate
from fenics import *

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


# ------------------------------ Trajectory computation -------------------- #

def compute_trajectory(integration_parameters_dict, particule_dict, mesh_dict, f, segments_list, zone):
    tmax = integration_parameters_dict['tmax']
    dt = integration_parameters_dict['dt']
    t = integration_parameters_dict['t']

    m_part = particule_dict['m_part']
    q_part = particule_dict['q_part']
    x0 = particule_dict['x0']
    y0 = particule_dict['y0']
    vx0 = particule_dict['vx0']
    vy0 = particule_dict['vy0']
    num_hole = particule_dict['num_hole']
    
    l_mot = mesh_dict['l_mot']
    l_vacuum = mesh_dict['l_vacuum']
    h_grid = mesh_dict['l_1'] + mesh_dict['l_2'] + mesh_dict['delta_vert_12']
    nb_hole = mesh_dict['nb_hole']
    
    Y=np.array([x0, y0, vx0, vy0])
    liste_x=[x0]
    liste_y=[y0]
    liste_vx=[vx0]
    liste_vy=[vy0]
    listet=[0] #pas un linspace au cas où on arrête avant
    N_impact=0

    while t<=tmax:
        k1=np.array(f(Y,t))
        k2=np.array(f(Y+.5*dt*k1, t+.5*dt))
        k3=np.array(f(Y+.5*dt*k2, t+.5*dt))
        k4=np.array(f(Y+dt*k3, t+dt))
        Y_pot=Y+dt*(1/6*k1+1/3*k2+1/3*k3+1/6*k4)

        while zone.inside(Point(Y_pot[0],Y_pot[1]))==False:
            xinter,yinter,n=coord_impact(Y[0],Y[1],Y_pot[0],Y_pot[1])
            if n=='xm' and num_hole+np.sign(xinter)<=nb_hole and num_hole+np.sign(xinter)>=0: #normale selon x avec miroir
                num_trou+=np.sign(xinter)
                Y,Y_pot=np.array([-xinter, -yinter, 0, 0]),\
                        np.array([-xinter-(xinter-Y_pot[0]), Y_pot[1], Y_pot[2], Y_pot[3]])

            elif n=='x' or n=='xm': #normale selon x sans miroir (ou avec mais en bout de chaine)
                Y,Y_pot=np.array([xinter, yinter, 0, 0]),\
                        np.array([xinter+(xinter-Y_pot[0]), Y_pot[1], -Y_pot[2], Y_pot[3]])
                N_impact+=1

            else:
                Y,Y_pot=np.array([xinter,yinter,0,0] ),\
                        np.array([Y_pot[0], yinter+(yinter-Y_pot[1]), Y_pot[2], -Y_pot[3]])
                N_impact+=1
            
        Y=Y_pot
        liste_x.append(Y[0])
        liste_y.append(Y[1])
        liste_vx.append(Y[2])
        liste_vy.append(Y[3])
        listet.append(t)
        t+=dt
        if Y[1]<-l_mot/2-h_grid-l_vacuum/2:
            break

    return liste_x, liste_y, liste_vx, liste_vy, N_impact
    

    

def intersection(x1,y1,x2,y2,x3,y3,x4,y4): #donne les coords d'inters des segments 1,2 et 3,4 #n='NO' si n'existe pas, 'x'ou'y' selon la normale #d'impact de  1,2 sur 3,4
    """
    Donne les coordonnées d'intersection des segments 1,2 et 3,4 en assumant que 3,4 est forcément vertical ou horizontal
    Renvoie aussi la normale au segment 3,4 (ici la normale au plan percuté)
    
    Si impact il y a eu, la sortie est de type (xi,yi,'x') ou (xi,yi,'y'). Sinon, renvoie (0,0,'NO')
    """
    if x3==x4: #segment vert
        if x2==x1: #parallelisme
            return (0,0,'NO')
        if (x1-x3)*(x2-x3)>0:  #du meme cote de ce segment
            return (0,0,'NO')
        pente_traj=(y2-y1)/(x2-x1)
        b_traj=y1-pente_traj*x1
        y_test=pente_traj*x3+b_traj
        if y_test<=max(y3,y4) and y_test>=min(y3,y4):
            return (x3,y_test,'x')
        return (0,0,'NO')
    else:    #segment horiz
        if y1==y2: #parallelisme
            return (0,0,'NO')
        if (y1-y3)*(y2-y3)>0:  #du meme cote de ce segment
            return (0,0,'NO')
        if x1==x2: #autre seg ss notion de pente
            if x1<=max(x3,x4) and x1>=min(x3,x4):
                return (x1,y3,'y')
            return (0,0,'NO')
        pente_traj=(y2-y1)/(x2-x1)
        b_traj=y1-pente_traj*x1
        x_test=(y3-b_traj)/pente_traj
        if x_test<=max(x3,x4) and x_test>=min(x3,x4):
            return (x_test,y3,'y')
        return (0,0,'NO')

def coord_impact(x1,y1,x2,y2, segments_list): #cherche où et selon quelle n le segments 1,2 coupe un bord
    """
    On sait que ce segment coupe un bord, On cherche alors quel segment il coupe ainsi que les coordonnées d'intersection et la normale
    Renvoie (xi,yi,n) ou n='x'ou'y'
    Par défaut (normalement jamais), renvoie (0,0,'NO') s'il ne coupe rien
    """
    for i in range(len(segments_list)):
        x3,y3,x4,y4=segments_list[i]
        xi,yi,n=intersection(x1,y1,x2,y2,x3,y3,x4,y4)
        if n!='NO':
            if i==0 or i==2 or i==6 or i==16:
                return (xi,yi,'xm')
            return (xi,yi,n)
    return (0,0,'NO')



# ------------------------- Bin ------------------- # 


"""
RK4 où Y est le vecteur au temps t de la boucle.
Les positions et vitesses et le temps correspondants sont en parallèles stockés dans les listes: liste_x, liste_vy ...
Notons une erreur notable, lors des calculs des Ki, on considère qu'il n'y a pas de rebonds possibles sur la trajectoire interpolée
d'où les cas de champ E nul dans f

Autre truc à étudier: si après le calcul du rebond ou du miroir, on retombe hors du mesh, potentiellement pas d'intersection à l'iteraion suivante
"""

'''while t<=tmax:
    k1=np.array(f(Y,t))
    k2=np.array(f(Y+.5*dt*k1, t+.5*dt))
    k3=np.array(f(Y+.5*dt*k2, t+.5*dt))
    k4=np.array(f(Y+dt*k3, t+dt))
    Y_pot=Y+dt*(1/6*k1+1/3*k2+1/3*k3+1/6*k4)
    if zone.inside(Point(Y_pot[0],Y_pot[1]))==True:
        Y=Y_pot
    else:
        xinter,yinter,n=coord_impact(Y[0],Y[1],Y_pot[0],Y_pot[1])
        N_impact+=1
        if n=='xm': #normale selon x avec miroir
            Y=np.array([-xinter-(xinter-Y_pot[0]),Y_pot[1],Y_pot[2],Y_pot[3]])
        if n=='x': #normale selon x
            Y=np.array([xinter+(xinter-Y_pot[0]),Y_pot[1],-Y_pot[2],Y_pot[3]])
        else:
            Y=np.array([Y_pot[0],yinter+(yinter-Y_pot[1]),Y_pot[2],-Y_pot[3]])
    liste_x.append(Y[0])
    liste_y.append(Y[1])
    liste_vx.append(Y[2])
    liste_vy.append(Y[3])
    listet.append(t)
    t+=dt
    if Y[1]<-l_mot/2-h_grid-l_vacuum/2:
        break'''
    
'''while t<=tmax:
    k1=np.array(f(Y,t))
    k2=np.array(f(Y+.5*dt*k1, t+.5*dt))
    k3=np.array(f(Y+.5*dt*k2, t+.5*dt))
    k4=np.array(f(Y+dt*k3, t+dt))
    Y_pot=Y+dt*(1/6*k1+1/3*k2+1/3*k3+1/6*k4)
    while zone.inside(Point(Y_pot[0],Y_pot[1]))==False:
        xinter,yinter,n=coord_impact(Y[0],Y[1],Y_pot[0],Y_pot[1], segments_list)
        if n=='xm': #normale selon x avec miroir
            Y,Y_pot=np.array([-xinter,-yinter,0,0]),np.array([-xinter-(xinter-Y_pot[0]),Y_pot[1],Y_pot[2],Y_pot[3]])
        elif n=='x': #normale selon x
            Y,Y_pot=np.array([xinter,yinter,0,0]),np.array([xinter+(xinter-Y_pot[0]),Y_pot[1],-Y_pot[2],Y_pot[3]])
            N_impact+=1
        else:
            Y,Y_pot=np.array([xinter,yinter,0,0]),np.array([Y_pot[0],yinter+(yinter-Y_pot[1]),Y_pot[2],-Y_pot[3]])
            N_impact+=1
    Y=Y_pot
    liste_x.append(Y[0])
    liste_y.append(Y[1])
    liste_vx.append(Y[2])
    liste_vy.append(Y[3])
    listet.append(t)
    t+=dt
    if Y[1]<-l_mot/2-h_grid-l_vacuum/2:
        break'''
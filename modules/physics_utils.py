from __future__ import print_function

# local import

from .particules import Particule
from .vector import MyVector


import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import dolfin
import mshr
import numpy as np
import scipy.integrate as integrate
from fenics import *

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
"""
class Particule:
    ""\"Classe définissant une particule caractérisée par :
    - son espèce
    - sa charge
    - sa masse
    - sa position x,y,z
    - sa vitesse vx, vy, vz
    ""\"
    def __init__(self, espece, q, m, x, y, z, vx, vy, vz):
        ""\"Constructeur""\"
        self.espece = espece
        self.q = q
        self.m = m
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
"""

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
        v=np.random.normal(3e3,1e3)
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

def One_step(liste_Y,n,segments_list,zone,mode_dict,mesh_dict,t,E,dt):
    Cond1=mode_dict['Elastique?']
    Cond2=mode_dict['Transfert de charge?']
    Cond3=mode_dict['Contact inter particules?']
    eta=mode_dict['perte u par contact'] #a priori idem pr I, I+ et I-
    p=mode_dict['proba perte q par contact'] #a priori idem pr I+ et I-

    if Cond3==True:
        print('error, Cond3 pas créée')
    
    particule=liste_Y[n][0]
    m=particule.get_mass()
    q=particule.get_charge()
    espece=particule.get_part_type()
    pos = particule.get_pos() # Vector object
    speed = particule.get_speed() 
    Y=np.array([pos.x, pos.y, pos.z, speed.x, speed.y, speed.z])
    k1=np.array(f(Y,t,m,q,zone,E))
    k2=np.array(f(Y+.5*dt*k1, t+.5*dt,m,q,zone,E))
    k3=np.array(f(Y+.5*dt*k2, t+.5*dt,m,q,zone,E))
    k4=np.array(f(Y+dt*k3, t+dt,m,q,zone,E))
    Y_pot=Y+dt*(1/6*k1+1/3*k2+1/3*k3+1/6*k4)

    while zone.inside(Point(Y_pot[0],Y_pot[1]))==False:
        xinter,yinter,n=coord_impact(Y[0],Y[1],Y_pot[0],Y_pot[1],segments_list)
        z=Y_pot[2]
        vz=Y_pot[5]
        
        if n=='x'or n=='xm': #normale selon x 
            #incidence = np.arctan((Y_pot[1]-Y[1])/(Y_pot[0]-Y[0]))
            #remplacer p par p*np.cos(incidence) et (1-eta) par (1-eta*np.cos(incidence)
            if Cond1==True:
                Y,Y_pot=np.array([xinter, yinter, z, 0, 0, vz]),\
                        np.array([xinter+(xinter-Y_pot[0]), Y_pot[1], z, -Y_pot[3], Y_pot[4], vz])
                if n=='x' and Cond2==True and q!=0 and np.random.random_sample()<=p: 
                    q,espece=0,'I'
            else:
                Y,Y_pot=np.array([xinter, yinter, z, 0, 0, vz]),\
                        np.array([xinter+(xinter-Y_pot[0]), Y_pot[1], z, -(1-eta)*Y_pot[3], (1-eta)*Y_pot[4], vz])
                if n=='x' and Cond2==True and q!=0 and np.random.random_sample()<=p:
                    q,espece=0,'I'
                    
        else: #normale selon y 
            #incidence = np.arctan((Y_pot[0]-Y[0])/(Y_pot[1]-Y[1]))
            #remplacer p par p*np.cos(incidence) et (1-eta) par (1-eta*np.cos(incidence)
            if Cond1==True:    
                Y,Y_pot=np.array([xinter,yinter,z, 0,0,vz] ),\
                        np.array([Y_pot[0], yinter+(yinter-Y_pot[1]), z, Y_pot[3], -Y_pot[4],vz])
                if Cond2==True and q!=0 and np.random.random_sample()<=p:
                    q,espece=0,'I'
            else:
                Y,Y_pot=np.array([xinter,yinter,z,0,0,vz] ),\
                        np.array([Y_pot[0], yinter+(yinter-Y_pot[1]),z, (1-eta)*Y_pot[3], -(1-eta)*Y_pot[4],vz])
                if Cond2==True and q!=0 and np.random.random_sample()<=p:
                    q,espece=0,'I'
    
    # update particule - what changes : charge, espece, 
    # possiblity mass even if it seems ignored there, position and speed.
    particule.set_pos(MyVector(Y_pot[0],Y_pot[1],Y_pot[2]))
    particule.set_speed(MyVector(Y_pot[3],Y_pot[4],Y_pot[5]))
    particule.set_charge(q) # here q is expected to be an int but that can be changed at any moment in class Particule
    particule.set_mass(m)
    particule.set_part_type(espece)
    
    return particule #Particule(espece, q, m, Y_pot[0], Y_pot[1], Y_pot[2], Y_pot[3], Y_pot[4], Y_pot[5])


# ------------------------------ Trajectory computation main functions -------------------- #

def compute_trajectory_bigN(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list, zone,E):
    """
    Renvoie la proportion d'espèce , 
    l'angle moyen, 
    et la norme de la vitesse moyenne en sortie.
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
    
    print("début d'initialisation")
    
    N2=int(p2*N)
    N3=int(p3*N)
    N1=N-N2-N3
    
    liste_Y=[]
    
    for n in range (N1):
        x,y,z,vx,vy,vz=distrib_init('I', mesh_dict)
        #liste_Y.append([Particule('I', 0, 127*u, x, y, z, vx, vy, vz),0])
        liste_Y.append([Particule(charge = 0, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I"),0]) # note that there still is a radius that can be set. By default it is set to IODINE_RADIUS.
    for n in range (N1,N1+N2):
        x,y,z,vx,vy,vz=distrib_init('I+', mesh_dict)
        #liste_Y.append([Particule('I+', e, 127*u, x, y, z, vx, vy, vz),0])
        liste_Y.append([Particule(charge = e, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I+"),0])
    for n in range (N1+N2,N):
        x,y,z,vx,vy,vz=distrib_init('I-', mesh_dict)
        #liste_Y.append([Particule('I-', -e, 127*u, x, y, z, vx, vy, vz),0])
        liste_Y.append([Particule(charge = -e, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I-"),0])
        
    np.random.shuffle(liste_Y)
    
    nombre_max_injecte_par_tour=int(DN*dt)+1 #au moins 1 si Dn trop petit
    nombre_tour_plein_debit=int(N/nombre_max_injecte_par_tour) 
    nombre_derniere_injection=N-nombre_tour_plein_debit*nombre_max_injecte_par_tour
    
    for i in range(nombre_max_injecte_par_tour, N):
        liste_Y[i][1]=-1
        
    ### On calcule les trajectoires en 2 temps
    ### On fait d'abord les tours où il y a injection au débit max en ne traitant qu'une partie des particules
    ### Puis on fait tourner en continu jusqu'à tmax ou à ce qu'il n'y ait plus de particules à l'intérieur
    
    print("début de l'injection")
        
    t=0
    Nb_out=0

    for i in range (nombre_tour_plein_debit):
        if np.random.random_sample()<max(dt/tmax,0.05):
            print ('Avancement de '+str(100*(t/tmax)%10)+'%, '+str(Nb_out)+' particules sorties')
        for n in range ((i+1)*nombre_max_injecte_par_tour):
            if liste_Y[n][1]!=1: 
                liste_Y[n][1]=0
                liste_Y[n][0]=One_step(liste_Y,n,segments_list,zone,mode_dict,mesh_dict,t,E,dt)
                if liste_Y[n][0].y<-l_mot/2-h_grid-l_vacuum/2:
                    liste_Y[n][1]=1
                    Nb_out+=1
        t+=dt
    
    print('toutes les particules sont injectées')
    
    while t<tmax and Nb_out<N:
        if np.random.random_sample()<max(dt/tmax,0.05):
            print ('Avancement de '+str(100*(t/tmax)%10)+'%, '+str(Nb_out)+' particules sorties')
        for n in range (N):
            if liste_Y[n][1]!=1: 
                liste_Y[n][1]=0
                liste_Y[n][0]=One_step(liste_Y,n,segments_list,zone,mode_dict,mesh_dict,t,E,dt)
                if liste_Y[n][0].get_pos().y<-l_mot/2-h_grid-l_vacuum/2:
                    liste_Y[n][1]=1
                    Nb_out+=1
        t+=dt
        
    if Nb_out==N:
        print("fin du calcul, critère d'arret: Nb_out")
    else:
        print("fin du calcul, critère d'arret: tmax")
    
    ### Traitement des données de sortie
    print('début traitement de données')
    
    N1f=0
    N2f=0
    N3f=0
    liste_vxf1=[]
    liste_vyf1=[]
    liste_vxf2=[]
    liste_vyf2=[]
    liste_vxf3=[]
    liste_vyf3=[]
   
    
    for n in range(N):
        if liste_Y[n][1]==1 and liste_Y[n][0].get_speed().y!=0:
            particule=liste_Y[n][0]
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
    liste_vxf1=np.array(liste_vxf1)
    liste_vyf1=np.array(liste_vyf1)
    liste_vxf2=np.array(liste_vxf2)
    liste_vyf2=np.array(liste_vyf2)
    liste_vxf3=np.array(liste_vxf3)
    liste_vyf3=np.array(liste_vyf3)
    liste_Vnorm1=np.power(np.power(liste_vxf1,2)+np.power(liste_vyf1,2), 1/2)
    liste_alpha1=-np.arctan(liste_vxf1/liste_vyf1)
    liste_Vnorm2=np.power(np.power(liste_vxf2,2)+np.power(liste_vyf2,2), 1/2)
    liste_alpha2=-np.arctan(liste_vxf2/liste_vyf2)
    liste_Vnorm3=np.power(np.power(liste_vxf3,2)+np.power(liste_vyf3,2), 1/2)
    liste_alpha3=-np.arctan(liste_vxf3/liste_vyf3)
    
    return [p1f,p2f,p3f], \
           [moyenne_amelioree(liste_alpha1),moyenne_amelioree(liste_alpha2),moyenne_amelioree(liste_alpha3)], \
           [moyenne_amelioree(liste_Vnorm1),moyenne_amelioree(liste_Vnorm2),moyenne_amelioree(liste_Vnorm3)]
    

def compute_trajectory_smallN(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list, zone,E):
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
    
    print("début d'initialisation")

    N2=int(p2*N)
    N3=int(p3*N)
    N1=N-N2-N3
    
    liste_Y=[]
    
    for n in range (N1):
        x,y,z,vx,vy,vz=distrib_init('I', mesh_dict)
        liste_Y.append([Particule(charge = 0, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I"),0]) # note that there still is a radius that can be set. By default it is set to IODINE_RADIUS.
    for n in range (N1,N1+N2):
        x,y,z,vx,vy,vz=distrib_init('I+', mesh_dict)
        liste_Y.append([Particule(charge = e, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I+"),0])
    for n in range (N1+N2,N):
        x,y,z,vx,vy,vz=distrib_init('I-', mesh_dict)
        liste_Y.append([Particule(charge = -e, mass = 127*u, pos = MyVector(x,y,z), speed = MyVector(vx,vy,vz), part_type = "I-"),0])
        
    np.random.shuffle(liste_Y)
 
    liste_t=[0]
    listes_x=[]
    listes_y=[]
    listes_vx=[]
    listes_vy=[]
    listes_q=[]
    for n in range(N):
        particule = liste_Y[n][0]
        pos = particule.get_pos()
        speed = particule.get_speed()
        listes_x.append([pos.x])
        listes_y.append([pos.y])
        listes_vx.append([speed.x])
        listes_vy.append([speed.y])
        listes_q.append([particule.get_charge()])

    nombre_max_injecte_par_tour=int(DN*dt)+1
    nombre_tour_plein_debit=int(N/nombre_max_injecte_par_tour)
    nombre_derniere_injection=N-nombre_tour_plein_debit*nombre_max_injecte_par_tour
    
    print("début de l'injection")
        
    for i in range(nombre_max_injecte_par_tour, N):
        liste_Y[i][1]=-1
        
    ### On calcule les trajectoires en 2 temps
    ### On fait d'abord les tours où il y a injection au débit max en ne traitant qu'une partie des particules
    ### Puis on fait tourner en continu jusqu'à tmax ou à ce qu'il n'y ait plus de particules à l'intérieur
        
    t=0
    Nb_out=0

    for i in range (nombre_tour_plein_debit):
        if np.random.random_sample()<max(dt/tmax,0.05):
            print ('Avancement de '+str(100*(t/tmax)%10)+'%, '+str(Nb_out)+' particules sorties')
        for n in range ((i+1)*nombre_max_injecte_par_tour):
            if liste_Y[n][1]!=1: 
                liste_Y[n][1]=0
                particule=One_step(liste_Y,n,segments_list,zone,mode_dict,mesh_dict,t,E,dt)
                liste_Y[n][0]=particule
                pos = particule.get_pos()
                speed = particule.get_speed()
                listes_x[n].append(pos.x)
                listes_y[n].append(pos.y)
                listes_vx[n].append(speed.x)
                listes_vy[n].append(speed.y)
                listes_q[n].append(particule.get_charge()) 
                if pos.y<-l_mot/2-h_grid-l_vacuum/2:
                    liste_Y[n][1]=1
                    Nb_out+=1
        t+=dt
        liste_t.append(t)
        
    print('toutes les particules sont injectées')
    
    
                 
    while t<tmax and Nb_out<N:
        if np.random.random_sample()<max(dt/tmax,0.05):
            print ('Avancement de '+str(100*(t/tmax)%10)+'%, '+str(Nb_out)+' particules sorties')
        for n in range (N):
            if liste_Y[n][1]!=1: 
                liste_Y[n][1]=0
                particule=One_step(liste_Y,n,segments_list,zone,mode_dict,mesh_dict,t,E,dt)
                liste_Y[n][0]=particule
                pos = particule.get_pos()
                speed = particule.get_speed()
                
                listes_x[n].append(pos.x)
                listes_y[n].append(pos.y)
                listes_vx[n].append(speed.x)
                listes_vy[n].append(speed.y)
                listes_q[n].append(particule.get_charge())
                
                if pos.y<-l_mot/2-h_grid-l_vacuum/2:
                    liste_Y[n][1]=1
                    Nb_out+=1
        t+=dt
        liste_t.append(t)
        
    if Nb_out==N:
        print("fin du calcul, critère d'arret: Nb_out")
    else:
        print("fin du calcul, critère d'arret: tmax")
    
    ### Traitement des données de sortie
    print('début traitement de données')
    
    N1f=0
    N2f=0
    N3f=0
    liste_vxf1=[]
    liste_vyf1=[]
    liste_vxf2=[]
    liste_vyf2=[]
    liste_vxf3=[]
    liste_vyf3=[]
   
    
     
    for n in range(N):
        if liste_Y[n][1]==1 and liste_Y[n][0].get_speed().y!=0:
            particule=liste_Y[n][0]
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
    liste_vxf1=np.array(liste_vxf1)
    liste_vyf1=np.array(liste_vyf1)
    liste_vxf2=np.array(liste_vxf2)
    liste_vyf2=np.array(liste_vyf2)
    liste_vxf3=np.array(liste_vxf3)
    liste_vyf3=np.array(liste_vyf3)
    liste_Vnorm1=np.power(np.power(liste_vxf1,2)+np.power(liste_vyf1,2), 1/2)
    liste_alpha1=-np.arctan(liste_vxf1/liste_vyf1)
    liste_Vnorm2=np.power(np.power(liste_vxf2,2)+np.power(liste_vyf2,2), 1/2)
    liste_alpha2=-np.arctan(liste_vxf2/liste_vyf2)
    liste_Vnorm3=np.power(np.power(liste_vxf3,2)+np.power(liste_vyf3,2), 1/2)
    liste_alpha3=-np.arctan(liste_vxf3/liste_vyf3)
    
    return [p1f,p2f,p3f], \
            [moyenne_amelioree(liste_alpha1),moyenne_amelioree(liste_alpha2),moyenne_amelioree(liste_alpha3)], \
            [moyenne_amelioree(liste_Vnorm1),moyenne_amelioree(liste_Vnorm2),moyenne_amelioree(liste_Vnorm3)], \
            listes_x, listes_y, listes_vx, listes_vy, listes_q, liste_t










    


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
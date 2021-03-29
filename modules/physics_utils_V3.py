from __future__ import print_function


import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import dolfin
import mshr
import numpy as np
from scipy.stats import maxwell
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
    L_3 = mesh_dict['L_3']
    l_3 = mesh_dict['l_3']
    delta_vert_12 = mesh_dict['delta_vert_12']
    delta_vert_23 = mesh_dict['delta_vert_23']
    mesh_resolution = mesh_dict['mesh_resolution']
    refine_mesh = mesh_dict['refine_mesh']

    h_grid = l_1 + l_2 + l_3 + delta_vert_12 + delta_vert_23

    class Boundary_top_mot(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] == l_mot/2
        
    class Boundary_bord_mot(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] < l_mot/2  and x[1] > - l_mot/2 + tol
        
    class Boundary_electrode1(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] <= - l_mot/2 + tol and x[1] >= - l_mot/2 - l_1 - tol
        
    class Boundary_inter_electrode(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            test_1 = on_boundary and x[1] < - l_mot/2 - l_1 - tol and x[1] > - l_mot/2 - l_1 - delta_vert_12 + tol
            test_2 = on_boundary and x[1] < - l_mot/2 - l_1 - delta_vert_12 - l_2 - tol and x[1] > - l_mot/2 - l_1 - delta_vert_12  - l_2 - delta_vert_23 + tol
            return test_1 or test_2
        
    class Boundary_electrode2(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] <= - l_mot/2 - l_1 - delta_vert_12 + tol and x[1] >= - l_mot/2 - \
                l_1 - delta_vert_12 - l_2 - tol and abs(x[0])<=L_mot/2
        
    class Boundary_electrode3(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14
            return on_boundary and x[1] <= - l_mot/2 - l_1 - delta_vert_12 - l_2 - delta_vert_23 + tol and x[1] >= - l_mot/2 - \
                l_1 - delta_vert_12 - l_2 - delta_vert_23 -l_3 - tol and abs(x[0])<=L_mot/2
        
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
    electrode3 = Boundary_electrode3()
    sup_vacuum = Boundary_sup_vacuum()
    inf_vacuum = Boundary_inf_vacuum()

    list_Phi = [phi_dict['Phi_top_mot'],phi_dict['Phi_bord_mot'],phi_dict['Phi_electrode1'], \
            phi_dict['Phi_inter_electrode'], phi_dict['Phi_electrode2'], phi_dict['Phi_electrode3'], \
            phi_dict['Phi_sup_vacuum'],phi_dict['Phi_inf_vacuum']]
    list_edges=[top_mot, bord_mot, electrode1, inter_electrode, electrode2, electrode3, sup_vacuum, inf_vacuum]
    bc=[]
    
    for i in range(len(list_edges)):
        if list_Phi[i]!='N':
            bc.append(DirichletBC(V, Constant(list_Phi[i]), list_edges[i]))


    u = TrialFunction(V)
    v = TestFunction(V)
    
    rhoelec=physics_consts_dict['rhoelec']
    l_rho=physics_consts_dict['l_rho']
    
    class Second_membre(UserExpression):
        def eval(self, value, x):
            ylim = -(0.5*l_mot + l_1 + delta_vert_12 + l_2 + delta_vert_23 + l_3)
            
            if l_rho!=0:
                
                if x[1] >= ylim + l_rho:
                    value[0] = rhoelec/physics_consts_dict['PERMITTIVITY']
                    
                elif x[1] >= ylim:
                    a = rhoelec*((-ylim - l_rho)**(-.5) - (-ylim)**(-.5))**(-1)
                    b = - a*(-ylim)**(-.5)
                    value[0] =  a * (-x[1])**(-.5) + b 
                    value[0] /= physics_consts_dict['PERMITTIVITY']
                    
                else:
                    value[0]=0
                    
            else:
                value[0]=rhoelec/physics_consts_dict['PERMITTIVITY']
    
    f = Second_membre(degree=2)

    #f = Constant(physics_consts_dict['rhoelec']/physics_consts_dict['PERMITTIVITY']) 
    a = dot(grad(u), grad(v))*dx
    L = f*v*dx
    Phi = Function(V)

    solve(a == L, Phi, bc)

    E = project(-grad(Phi), W)

    return Phi,E,f


# ------------------------------ Trajectory computation auxiliary functions -------------------- #

class Particule:
    """Classe définissant une particule caractérisée par :
    - son espèce
    - sa charge
    - sa masse
    - sa position x,y,z
    - sa vitesse vx, vy, vz
    """
    def __init__(self, espece, q, m, x, y, z, vx, vy, vz):
        """Constructeur"""
        self.espece = espece
        self.q = q
        self.m = m
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        
def moyenne_amelioree(liste):
    if len(liste)==0:
        return 0
    return np.mean(liste)

def ecart_type_ameliore(liste):
    if len(liste)==0:
        return 0
    return np.std(liste)

def distrib_init(mesh_dict, injection_dict, mode_dict,physics_consts_dict):
    """
    Renvoie x,y,z aléatoirement en suivant une distribution propre à l'espèce
    """
    L_mot = mesh_dict['L_mot']
    l_mot = mesh_dict['l_mot']
    L_1 = mesh_dict['L_1']
    l_1 = mesh_dict['l_1']
    gamma=injection_dict['gamma']
    Tion=injection_dict['Tion']
    cond=mode_dict['X0 fixe ?']
    u=physics_consts_dict['M_NUCLEON']
    Kb=physics_consts_dict['K_BOLTZMANN']
    tol=1e-5
    loc_ion=0
    scale_ion=np.sqrt(Kb*Tion/127*u)
    
    alpha_lim=np.arctan((L_mot-L_1)/(2*gamma*l_1))
    if cond==False:
        alpha=np.random.uniform(np.pi/2-alpha_lim+tol,np.pi/2+alpha_lim-tol)
        v=maxwell.rvs(loc=loc_ion, scale=scale_ion)
    else:
        alpha=np.pi/2-alpha_lim/3
        v=np.sqrt(8*Kb*Tion/(np.pi*127*u))
    R=np.sqrt((gamma*l_1)**2+((L_mot-L_1)/2)**2)
    x=R*np.cos(alpha)
    y=-l_mot/2-gamma*l_1+R*np.sin(alpha)
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

def coord_impact(x1,y1,x2,y2, segments_list, last_impact): #cherche où et selon quelle n le segments 1,2 coupe un bord
    """
    On sait que ce segment coupe un bord, On cherche alors quel segment il coupe ainsi que les coordonnées d'intersection, la normale et le numéro du segment en question. L'argument last_impact est éventuellement le dernier segment impacté, qu'il ne faut pas reprendre en compte
    Renvoie (xi,yi,n,i) ou n='x'ou'y'
    Par défaut (normalement jamais), renvoie (0,0,'NO',None) s'il ne coupe rien
    """
    for i in range(len(segments_list)):
        if i!=last_impact :
            x3,y3,x4,y4=segments_list[i]
            xi,yi,n=intersection(x1,y1,x2,y2,x3,y3,x4,y4)
            if n!='NO':
                if i==0 or i==2 or i==6 or i==10 or i==20 or i==24:
                    return (xi,yi,'xm',i)
                return (xi,yi,n,i)
    return (0,0,'NO',None)

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

def One_step(liste_Y,n,segments_list,zone,mode_dict,mesh_dict,physics_consts_dict,t,E,dt, injection_dict):
    """
    liste_Y est la liste des couples [particule,i] de la simulation, n est la particule à faire avancer.
    renvoie la particule dans son nouvel état.
    """
    Cond1=mode_dict['Choc symétrique sur parois ?']
    Cond2=mode_dict['Contact inter particules ?']
    Cond3=mode_dict['Perte de charge sur parois ?']
    Cond4=mode_dict['Perte d energie sur parois ?']
    d1,d2=physics_consts_dict['param densité neutre']
    K0_bis,K1_bis,K2_bis = mode_dict['param moy angle rebond']
    A_bis,B_bis = mode_dict['param sigma angle rebond']
    K0,K1,K2 = mode_dict['param moy energie rebond']
    A,B = mode_dict['param sigma energie rebond']
    e=physics_consts_dict['CHARGE']
    u=physics_consts_dict['M_NUCLEON']
    Kb=physics_consts_dict['K_BOLTZMANN']
    l_mot=mesh_dict['l_mot']
    Tneutre=injection_dict['Tneutre']
    
    loc_neutre=0
    scale_neutre=np.sqrt(Kb*Tneutre/127*u)
    
    particule=liste_Y[n][0]
    m=particule.m
    q=particule.q
    espece=particule.espece
    
    Y=np.array([particule.x, particule.y, particule.z, particule.vx, particule.vy, particule.vz])
    k1=np.array(f(Y,t,m,q,zone,E))
    k2=np.array(f(Y+.5*dt*k1, t+.5*dt,m,q,zone,E))
    k3=np.array(f(Y+.5*dt*k2, t+.5*dt,m,q,zone,E))
    k4=np.array(f(Y+dt*k3, t+dt,m,q,zone,E))
    Y_pot=Y+dt*(1/6*k1+1/3*k2+1/3*k3+1/6*k4)
    
    last_impact=None

    while zone.inside(Point(Y_pot[0],Y_pot[1]))==False:
        xinter,yinter,n,last_impact=coord_impact(Y[0],Y[1],Y_pot[0],Y_pot[1],segments_list,last_impact)
        z=Y_pot[2]
        vz=Y_pot[5]
        
##############################################################################################################################
        if n=='xm': #normale selon x sans choc
            Y,Y_pot=np.array([xinter, 
                              yinter, 
                              z, 
                              0, 
                              0, 
                              vz]),\
                    np.array([xinter+(xinter-Y_pot[0]), 
                              Y_pot[1], 
                              z, 
                              -Y_pot[3], 
                              Y_pot[4], 
                              vz])
##############################################################################################################################
        elif n=='x': #normale selon x avec choc
        
            sigmax=np.sign(Y_pot[0]-Y[0])
            sigmay=np.sign(Y_pot[1]-Y[1])
            incidence = sigmax*sigmay*np.arctan((Y_pot[1]-Y[1])/(Y_pot[0]-Y[0]))
            moy_E=K0+K1*np.exp(180*incidence/(np.pi*K2))
            sigma_E=A*180*incidence/np.pi + B
            coef_perte_energie=1-moy_E+np.random.random_sample()*sigma_E
            coef_perte_energie*=int(Cond4)
            
            if Cond1==True : #choc symetrique
                Y,Y_pot=np.array([xinter, 
                                  yinter, 
                                  z, 
                                  0, 
                                  0, 
                                  vz]),\
                        np.array([xinter+(xinter-Y_pot[0]), 
                                  Y_pot[1], 
                                  z, 
                                  -np.sqrt(1-coef_perte_energie)*Y_pot[3], 
                                  np.sqrt(1-coef_perte_energie)*Y_pot[4], 
                                  np.sqrt(1-coef_perte_energie)*vz])
                
            else: #choc asymetrique
                Vit_xy=np.sqrt(1-coef_perte_energie)*np.sqrt(Y_pot[3]**2+Y_pot[4]**2)
                Pos_xy=np.sqrt((xinter-Y_pot[0])**2+(Y_pot[1]-yinter)**2)
                moy_angle=K0_bis+K1_bis*np.exp(180*incidence/(np.pi*K2_bis))
                sigma_angle=A_bis*180*incidence/np.pi + B_bis
                i_prime=moy_angle+np.random.random_sample()*sigma_angle
                i_prime*=np.pi/180
                Y,Y_pot=np.array([xinter, 
                                  yinter, 
                                  z, 
                                  0, 
                                  0, 
                                  vz]),\
                        np.array([xinter-sigmax*np.cos(i_prime)*Pos_xy, 
                                  yinter+sigmay*np.sin(i_prime)*Pos_xy, 
                                  z, 
                                  -sigmax*np.cos(i_prime)*Vit_xy, 
                                  sigmay*np.sin(i_prime)*Vit_xy, 
                                  np.sqrt(1-coef_perte_energie)*vz])
                
            if q!=0 and np.random.random_sample()<=np.cos(incidence) and Cond3==True:
                q,espece=0,'I'
                
##############################################################################################################################
        else: #normale selon y 
        
            sigmax=np.sign(Y_pot[0]-Y[0])
            sigmay=np.sign(Y_pot[1]-Y[1])
            incidence = sigmax*sigmay*np.arctan((Y_pot[0]-Y[0])/(Y_pot[1]-Y[1]))
            moy_E=K0+K1*np.exp(180*incidence/(np.pi*K2))
            sigma_E=A*180*incidence/np.pi + B
            coef_perte_energie=1-moy_E+np.random.random_sample()*sigma_E
            coef_perte_energie*=int(Cond4)
            
            if Cond1==True: #choc symetrique
                Y,Y_pot=np.array([xinter,
                                  yinter,
                                  z,
                                  0,
                                  0,
                                  vz] ),\
                        np.array([Y_pot[0], 
                                  yinter+(yinter-Y_pot[1]),
                                  z, 
                                  np.sqrt(1-coef_perte_energie)*Y_pot[3], 
                                  -np.sqrt(1-coef_perte_energie)*Y_pot[4],
                                  np.sqrt(1-coef_perte_energie)*vz])
                
            else: #choc asymetrique
                Vit_xy=np.sqrt(1-coef_perte_energie)*np.sqrt(Y_pot[3]**2+Y_pot[4]**2)
                Pos_xy=np.sqrt((Y_pot[0]-xinter)**2+(yinter-Y_pot[1])**2)
                moy_angle=K0_bis+K1_bis*np.exp(180*incidence/(np.pi*K2_bis))
                sigma_angle=A_bis*180*incidence/np.pi + B_bis
                i_prime=moy_angle+np.random.random_sample()*sigma_angle
                i_prime*=np.pi/180
                Y,Y_pot=np.array([xinter,
                                  yinter,
                                  z,
                                  0,
                                  0,
                                  vz] ),\
                        np.array([xinter+sigmax*np.sin(i_prime)*Pos_xy, 
                                  yinter-sigmay*np.cos(i_prime)*Pos_xy, 
                                  z, 
                                  sigmax*np.sin(i_prime)*Vit_xy, 
                                  -sigmay*np.cos(i_prime)*Vit_xy, 
                                  np.sqrt(1-coef_perte_energie)*vz])
                
            if q!=0 and np.random.random_sample()<=np.cos(incidence) and Cond3==True:
                q,espece=0,'I'
                
    sect_eff=1e-19
##############################################################################################################################
    V=np.sqrt(Y_pot[3]**2+Y_pot[4]**2+Y_pot[5]**2)
    if Y_pot[1]<-l_mot/2:
        d=d1
    else:
        d=d2
    if np.random.random_sample()<dt*d*sect_eff*V and Cond2==True and espece=='I-':
        q,espece=0,'I'
        vx_ini=np.random.normal(loc_neutre,scale_neutre)
        vy_ini=np.random.normal(loc_neutre,scale_neutre)
        vz_ini=np.random.normal(loc_neutre,scale_neutre)
        liste_Y.append([Particule('I-', -e, 127*u, Y_pot[0], Y_pot[1], Y_pot[2], vx_ini, vy_ini, vz_ini),.5])
        
    return Particule(espece, q, m, Y_pot[0], Y_pot[1], Y_pot[2], Y_pot[3], Y_pot[4], Y_pot[5]),liste_Y

# ------------------------------ Trajectory computation main function -------------------- #

def compute_trajectory(integration_parameters_dict, injection_dict, mesh_dict, mode_dict, segments_list, zone, E, physics_consts_dict, details, inside):
    """
    Renvoie la proportion d'espèce , 
    l'angle moyen, 
    et la norme de la vitesse moyenne en sortie.
    
    Ainsi que la liste des instants t et les listes de liste de x,y,vx,vy,q de chaque particule
    """
    
    tmax=integration_parameters_dict['tmax']
    dt=integration_parameters_dict['dt']
    percent=integration_parameters_dict['%Nout']
    
    e=physics_consts_dict['CHARGE']
    u=physics_consts_dict['M_NUCLEON']
    
    l_mot=mesh_dict['l_mot']
    l_1=mesh_dict['l_1']
    l_2=mesh_dict['l_2']
    l_3=mesh_dict['l_3']
    delta_vert_12=mesh_dict['delta_vert_12']
    delta_vert_23=mesh_dict['delta_vert_23']
    l_vacuum=mesh_dict['l_vacuum']
    h_grid=l_1 + l_2 + l_3 + delta_vert_12 + delta_vert_23
    
    N=injection_dict['Nombre de particules I- en entrée']
    DN=injection_dict['débit de particule en entrée']

##############################################################################################################################
    ### On initialise une liste_Y comportant N lignes, chaque ligne correspond à une particule et est un couple (Particule,Entier)
    ### l'entier vaut -1 si la particule n'est pas encore injectée, 0 si est dedans et 1 si est sortie
    print("début d'initialisation")
    liste_Y=[]
    
    for n in range (N):
        x,y,z,vx,vy,vz=distrib_init(mesh_dict,injection_dict, mode_dict,physics_consts_dict)
        liste_Y.append([Particule('I-', -e, 127*u, x, y, z, vx, vy, vz),0])
   
    if details==True:
        liste_t=[0]
        listes_x=[]
        listes_y=[]
        listes_vx=[]
        listes_vy=[]
        listes_vz=[]
        listes_q=[]
        for n in range(N):
            listes_x.append([liste_Y[n][0].x])
            listes_y.append([liste_Y[n][0].y])
            listes_vx.append([liste_Y[n][0].vx])
            listes_vy.append([liste_Y[n][0].vy])
            listes_vz.append([liste_Y[n][0].vz])
            listes_q.append([liste_Y[n][0].q])
      
 ############################################################################################################################## 
    ### On calcule les trajectoires en 2 temps
    ### On fait d'abord les tours où il y a injection au débit max en ne traitant qu'une partie des particules
    ### Puis on fait tourner en continu jusqu'à tmax ou à ce qu'il n'y ait plus de particules à l'intérieur
    print("début de l'injection")
        
    t=0
    Nb_out=0
    last_percent=-1
    nb_nv_ion=0
                 
    while t<tmax and Nb_out<int(percent*N/100):
        if round(100*(t/tmax))!=last_percent:
            print ('Avancement de '+str(round(100*(t/tmax)))+'%, '+str(Nb_out)+' particules sorties, '+str(nb_nv_ion)+' nouveaux ions')
            last_percent=round(100*(t/tmax))
        for n in range (len(liste_Y)):
            #print('calcul de la part num'+str(n))
            if liste_Y[n][1]!=1:
                nb_avant=len(liste_Y)
                particule,liste_Y=One_step(liste_Y,n,segments_list,zone,mode_dict,mesh_dict,physics_consts_dict,t,E,dt, injection_dict)
                liste_Y[n][0]=particule
                if len(liste_Y)!=nb_avant and details==True:
                    nb_nv_ion+=1
                    new_particule=liste_Y[-1][0]
                    nb_etapes_passees=len(listes_x[n])
                    
                    listes_x.append([0]*nb_etapes_passees)
                    listes_y.append([0]*nb_etapes_passees)
                    listes_vx.append([0]*nb_etapes_passees)
                    listes_vy.append([0]*nb_etapes_passees)
                    listes_vz.append([0]*nb_etapes_passees)
                    listes_q.append([0]*nb_etapes_passees)
                    listes_x[-1].append(new_particule.x)
                    listes_y[-1].append(new_particule.y)
                    listes_vx[-1].append(new_particule.vx)
                    listes_vy[-1].append(new_particule.vy)
                    listes_vz[-1].append(new_particule.vz)
                    listes_q[-1].append(new_particule.q)
                    
                if details==True:
                    listes_x[n].append(particule.x)
                    listes_y[n].append(particule.y)
                    listes_vx[n].append(particule.vx)
                    listes_vy[n].append(particule.vy)
                    listes_vz[n].append(particule.vz)
                    listes_q[n].append(particule.q)
                
                if particule.y<-l_mot/2-h_grid-l_vacuum/2:
                    liste_Y[n][1]=1
                    Nb_out+=1
        t+=dt
        if details==True:
            liste_t.append(t)
            
##############################################################################################################################
    if Nb_out>=int(percent*N/100):
        print("fin du calcul, critère d'arret: % Nb_out, t="+str(t))
    else:
        print("fin du calcul, critère d'arret: tmax, N_out="+str(Nb_out))
    
    ### Traitement des données de sortie
    print('début traitement de données')
    
    Nf_ion=0
    Nf_neutre=0
    vxf_ion=[]
    vyf_ion=[]
    vzf_ion=[]
    vxf_neutre=[]
    vyf_neutre=[]
    vzf_neutre=[]
    
    if inside==True:
        vxf_ion_inside=[]
        vyf_ion_inside=[]
        vzf_ion_inside=[]
        vxf_neutre_inside=[]
        vyf_neutre_inside=[]
        vzf_neutre_inside=[]
   
    for n in range(N):
        if liste_Y[n][1]==1 and liste_Y[n][0].vy!=0: 
            particule=liste_Y[n][0]
            if particule.espece=='I':
                Nf_neutre+=1
                vxf_neutre.append(particule.vx)
                vyf_neutre.append(particule.vy)
                vzf_neutre.append(particule.vz)
            else:
                Nf_ion+=1
                vxf_ion.append(particule.vx)
                vyf_ion.append(particule.vy)
                vzf_ion.append(particule.vz)
        elif liste_Y[n][1]==0 and inside==True:
            particule=liste_Y[n][0]
            if particule.espece=='I':
                vxf_neutre_inside.append(particule.vx)
                vyf_neutre_inside.append(particule.vy)
                vzf_neutre_inside.append(particule.vz)
            else:
                vxf_ion_inside.append(particule.vx)
                vyf_ion_inside.append(particule.vy)
                vzf_ion_inside.append(particule.vz)
    
    if Nb_out!=0:
        pf_ion=Nf_ion/Nb_out
        pf_neutre=Nf_neutre/Nb_out
    else:
        pf_ion=0
        pf_neutre=0
        
    vxf_ion=np.array(vxf_ion)
    vyf_ion=np.array(vyf_ion)
    vzf_ion=np.array(vzf_ion)
    vxf_neutre=np.array(vxf_neutre)
    vyf_neutre=np.array(vyf_neutre)
    vzf_neutre=np.array(vzf_neutre)
    
    if inside==True:
        vxf_ion_inside=np.array(vxf_ion_inside)
        vyf_ion_inside=np.array(vyf_ion_inside)
        vzf_ion_inside=np.array(vzf_ion_inside)
        vxf_neutre_inside=np.array(vxf_neutre_inside)
        vyf_neutre_inside=np.array(vyf_neutre_inside)
        vzf_neutre_inside=np.array(vzf_neutre_inside)
    
    Vf_ion=np.power(np.power(vxf_ion,2)+np.power(vyf_ion,2)+np.power(vzf_ion,2), 1/2)
    anglef_ion=-np.arctan(vxf_ion/vyf_ion)
    Vf_neutre=np.power(np.power(vxf_neutre,2)+np.power(vyf_neutre,2)+np.power(vzf_neutre,2), 1/2)
    anglef_neutre=-np.arctan(vxf_neutre/vyf_neutre)
    
    if inside==True:
        Vf_ion_inside=np.power(np.power(vxf_ion_inside,2)+np.power(vyf_ion_inside,2)+np.power(vzf_ion_inside,2), 1/2)
        Vf_neutre_inside=np.power(np.power(vxf_neutre_inside,2)+np.power(vyf_neutre_inside,2)+np.power(vzf_neutre_inside,2), 1/2)
    
    if details==True:
        if inside==True:
            return Nb_out,[pf_neutre,pf_ion], \
                    [moyenne_amelioree(anglef_neutre),moyenne_amelioree(anglef_ion)], \
                    [ecart_type_ameliore(anglef_neutre),ecart_type_ameliore(anglef_ion)], \
                    [moyenne_amelioree(Vf_neutre),moyenne_amelioree(Vf_ion)], \
                    [ecart_type_ameliore(Vf_neutre),ecart_type_ameliore(Vf_ion)], \
                    listes_x, listes_y, listes_vx, listes_vy, listes_q, liste_t,\
                    anglef_neutre, anglef_ion, Vf_neutre, Vf_ion,\
                    Vf_neutre_inside, Vf_ion_inside
        else:
            return Nb_out,[pf_neutre,pf_ion], \
                    [moyenne_amelioree(anglef_neutre),moyenne_amelioree(anglef_ion)], \
                    [ecart_type_ameliore(anglef_neutre),ecart_type_ameliore(anglef_ion)], \
                    [moyenne_amelioree(Vf_neutre),moyenne_amelioree(Vf_ion)], \
                    [ecart_type_ameliore(Vf_neutre),ecart_type_ameliore(Vf_ion)], \
                    listes_x, listes_y, listes_vx, listes_vy, listes_q, liste_t,\
                    anglef_neutre, anglef_ion, Vf_neutre, Vf_ion
    else:
        if inside==True:
            print('erreur, details=True and inside=False impossible')
        else:
            return Nb_out,[pf_neutre,pf_ion], \
                    [moyenne_amelioree(anglef_neutre),moyenne_amelioree(anglef_ion)], \
                    [ecart_type_ameliore(anglef_neutre),ecart_type_ameliore(anglef_ion)], \
                    [moyenne_amelioree(Vf_neutre),moyenne_amelioree(Vf_ion)], \
                    [ecart_type_ameliore(Vf_neutre),ecart_type_ameliore(Vf_ion)]

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
from modules.particles_init import get_particles
#from modules.plotting_utils import ?

mesh_dict = {
    'L_mot' : .01,
    'l_mot' : .003,
    'L_1' : .0045,
    'l_1': .005,
    'L_2' : .007,
    'l_2' : .015,
    'delta_vert_12' : .005,
    'L_vacuum' : .01, # x
    'l_vacuum': .005, # y
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

# ------------------------ My part ------------------------ #
# my modules imports
from modules.collisions_handler import CollisionHandler
from modules.particules import Particule
from modules.utils import segment
from modules.vector import MyVector
from modules.grid import Grid
from modules.utils import get_mass_part
from modules.integration_schemes import scipy_integrate_solve_ivp, rk4, euler_explicit

# import for analysis
from modules.data_analysis import DataSaver, DataAnalyser

# imports
from dolfin import Point
from random import random
import numpy as np
from scipy.stats import maxwell, norm

# plotting
from tqdm import tqdm

# comparaison
from time import time
# ----------------- saving test to file ? --------------- #
saving_directory = 'tests/' # be careful, it won't check if this the directory is already created ...
tests_summary_file_name = "tests_summary_complete_1"
save_test = False
saving_period = 10
# ----------------- id simulation -------------------- #
id_test = 1
# --------------- Default analysis ? ---------------- #
perform_default_analysis = False
#----------------- debug parameters --------------------#
debug = True
test_config = True

#----------------- physics properties --------------------#

real_particle_density = 1e20 # for I
effective_diameter = 4e-10 # roughly the diameter of I
max_speed = 30000 # 20000 m/s we say ~> Maxwellian distribution
# TODO investigate max_speed importance

factor_in_front_of_mean_free_path = 10 # 1 mean free path = factor_in_front_of_mean_free_path cells
mean_free_path = 1/(np.sqrt(2)*np.pi*effective_diameter*effective_diameter*real_particle_density)
mean_free_time = mean_free_path/max_speed
dt = 0.25 * mean_free_time / factor_in_front_of_mean_free_path

mean_particles_number_per_cell = 50

MAX_INTEGRATION_STEP = 300
#----------------- Space properties --------------------#

# TODO we need to use the space
# resolution along each previous dimension

# rectangle of size l1*l2
segment_X_list = []
segment_Y_list = []
for segment in segments_list:
    segment_X_list.append(segment[0])
    segment_X_list.append(segment[2])
    segment_Y_list.append(segment[1])
    segment_Y_list.append(segment[3])
max_x, min_x = max(segment_X_list), min(segment_X_list)
max_y, min_y = max(segment_Y_list), min(segment_Y_list)

l1,l2 = max_x - min_x, max_y - min_y

res1, res2 = int(factor_in_front_of_mean_free_path*l1/mean_free_path), int(factor_in_front_of_mean_free_path*l2/mean_free_path)
nb_cells = res1*res2
if(debug): print("Dimension : {}x{} m".format(round(l1,2),round(l2,2)))

l3 = 0.1 # 10 cm ?

#----------------- Grid creation ----------------------#
dtype = "LinkedList"

# TODO : make sure it works
offset_x = -min_x
offset_y = -min_y
my_grid = Grid(l1,l2,[res1,res2], offsets = [offset_x,offset_y], dtype = dtype) # LinkedList , DynamicArray

# should change how the grid works...

#--------------- Particles creation -------------------#
N_particles_real = int(real_particle_density*l1*l2*l3) # this is the REAL number of particles
N_particles_simu = int(mean_particles_number_per_cell*res1*res2)
Ne = int(N_particles_real/N_particles_simu)

if(debug): print("There is {} particles in the simulation. One particle accounts for {} real particles.".format("{:e}".format(N_particles_simu),"{:e}".format(Ne)))

types = ['I']
numbers = [N_particles_simu]
speed_init_type = ['maxwellian']
speed_init_params = [[500,3000]]
list_particles = get_particles(types, numbers, speed_init_type, speed_init_params, effective_diameter, zone, [offset_x, offset_y], [l1,l2], verbose = False, debug = False)

#--------------- Rectangle creation -------------------#

# this should be changed completely !
rectangle_walls = segments_list

#[segment(Point(0,0),Point(0,l2)), segment(Point(0,0),Point(l1,0)), \
#    segment(Point(l1,0),Point(l1,l2)), segment(Point(0,l2),Point(l1,l2))]

#---------- Adding particles to the grid -------------#
for particle in list_particles:
    my_grid.add(particle)

#---------- Creating collision handler ---------------#

    # update function*
"""
def f(Y,t,m,q):
    vx=Y[3]
    vy=Y[4]
    vz=Y[5]
    
    ax = 0 
    ay = 0
    az = 0
    return np.array([vx, vy, vz, ax, ay, az])
"""

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
    return np.array([vx, vy, vz, ax, ay, az])

    # parameters
eta = 0
p = 0

    # DSMC params
DSMC_params = {
    'vr_max' : 2*max_speed,
    'effective_diameter':  effective_diameter,
    'Ne' : Ne, # this is the number of real particles one simulated particle represents.
    'cell_volume' : l1*l2*l3/(res1*res2) # since we are in 2D I don't really know what to add here actually... For now, I add the 3rd dimension rough size, that is l3
}
use_particles_collisions = False
use_DSMC = True
integration_scheme = euler_explicit # scipy_integrate_solve_ivp , rk4 , euler_explicit


# we now just have to place the 1 at the right place !
sparsed_space  = my_grid.fill_sparsed_space_from_initial_particles_position(list_particles)

if(debug):
    print("Size of the space : {}x{}, resolutions : {}x{}, offsets : {}x{}".format(round(l1,2),round(l2,2),res1,res2,offset_x,offset_y))

    print("Sparsed matrix : {}".format(sparsed_space))
    
    fig, ax = plt.subplots(figsize=(15,10))
    fig.suptitle("Working space")
    for k in range(len(segments_list)):
        X = [segment_X_list[2*k],segment_X_list[2*k+1]]
        Y = [segment_Y_list[2*k],segment_Y_list[2*k+1]]
        plt.plot(Y,X)
    plt.plot()
    plt.show()

collisionHandler = CollisionHandler(list_particles, rectangle_walls, f, eta, p, use_particles_collisions = use_particles_collisions, \
    use_DSMC = use_DSMC, grid = my_grid, DSMC_params = DSMC_params, sparsed_space = sparsed_space, integration_scheme=integration_scheme)

#------------------- Simulation -------------------#
if(save_test):
    params_dict = {
        'id_test' : id_test,
        'total_number_of_particles' : N_particles_simu,
        'path_to_data' : saving_directory+str(id_test),
        'real_particle_density' : real_particle_density,
        'effective_diameter' : effective_diameter,
        'max_speed' : max_speed,
        'mean_particles_number_per_cell' : mean_particles_number_per_cell,
        'nb_cells' : nb_cells,
        'mean_free_path' : mean_free_path,
        'mean_free_time' : mean_free_time,
        'dt' : dt,
        'MAX_INTEGRATION_STEP' : MAX_INTEGRATION_STEP,
        'nb_I' : N_particles_simu,
        'nb_I_real' : N_particles_real,
        'Ne' : Ne,
        'eta': eta,
        'loss_charge_proba' : p,
        'use_particles_collisions' : use_particles_collisions,
        'use_DSMC' : use_DSMC,
        'integration_scheme' : integration_scheme.__name__,
        'speed_init_type' : speed_init_type,
        'saving_period' : saving_period
    }
    data_analyser = DataSaver(list_particles, name_test = str(id_test), saving_directory = saving_directory)
    #data_analyser.save_test_params(tests_summary_file_name, params_dict, use_saving_directory = False)
    # integration params
t = 0

if(debug): print("\nSTARTING SIMULATION...\n")
# simulation

elapsed_time = time()

if(save_test):
        data_analyser.save_everything_to_one_csv()

if(not test_config):
    for k in tqdm(range(MAX_INTEGRATION_STEP)):
        collisionHandler.step(dt, t, [zone,E])
        t+=dt
        if(save_test and (k%saving_period==0 or k == MAX_INTEGRATION_STEP-1)): # we are saving last frame anyway
            data_analyser.save_everything_to_one_csv()

    if(debug): print("\nElapsed  time for {} iterations with {} particules and with {} data structure : {}".format(MAX_INTEGRATION_STEP, N_particles_simu, dtype, round(time()-elapsed_time,3)))


    if(debug): print("\nNumber of collisions : {}".format("{:e}".format(collisionHandler.get_collisions_count())))
    
    data_analyser.save_test_params(tests_summary_file_name, params_dict, use_saving_directory = False)
    
    if(perform_default_analysis):
        data_analyser = DataAnalyser(tests_summary_file_name)
        data_analyser.load_test(id_test)
        data_analyser.draw_particles()
        data_analyser.draw_hist_distribution('vx')
        data_analyser.draw_hist_distribution('vy')
        data_analyser.draw_hist_distribution('speed_norm_squared')
        data_analyser.draw_hist_distribution('speed_norm')
        data_analyser.draw_spatial_distribution(None, vmin = 0, vmax = 1.5*mean_particles_number_per_cell)
        data_analyser.draw_spatial_distribution('vx', vmin = -5e3, vmax = 5e3)
        data_analyser.draw_spatial_distribution('vy', vmin = -5e3, vmax = 5e3)
        data_analyser.draw_spatial_distribution('speed_norm_squared', vmin = 1e6, vmax = 25e6)
        data_analyser.draw_spatial_distribution('speed_norm', vmin = 1e3, vmax = 5e3)
else: 
    if(debug): print("Testing one step computation ...", end = " ")
    collisionHandler.step(dt, t, [zone,E])
    if(debug) : print("\t[OK]")
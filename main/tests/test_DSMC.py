from main import Grid, get_mass_part, DataAnalyser, \
     DataSaver, get_particles, segment, CollisionHandler, euler_explicit
from numpy.lib.function_base import delete
# imports
from dolfin import Point
from random import random
import numpy as np
from scipy.stats import maxwell, norm
import csv
from os import path, mkdir

# plotting
from tqdm import tqdm

# comparaison
from time import time
# ----------------- saving test to file ? --------------- #
saving_directory = 'C/' # be careful, it won't check if this the directory is already created ...
# checking if the file exists
if(not path.exists(saving_directory)):
    try:
        mkdir(saving_directory)
        print('Successfully created the saving directory {}'.format(saving_directory))
    except:
        print('Could not create the saving directory {}'.format(saving_directory))

tests_summary_file_name = "C"
save_test = True
saving_period = 1
# ----------------- id test -------------------- #
id_test = 1
# --------------- Default analysis ? ---------------- #
perform_default_analysis = False
#----------------- debug parameters --------------------#
debug = True
test_config = False 

#----------------- physics properties --------------------#

real_particle_density = 1e20 # for I
effective_diameter = 4e-10 # roughly the diameter of I
max_speed = 600 # 20000 m/s we say ~> Maxwellian distribution
mean_speed = 300
# TODO investigate max_speed importance

mean_free_path = 1/(np.sqrt(2)*np.pi*effective_diameter*effective_diameter*real_particle_density)
mean_free_time = mean_free_path/mean_speed
dt = 0.25 * mean_free_time # may be take a smaller one this one is kind of big still.

mean_particles_number_per_cell = 100
MAX_INTEGRATION_STEP = 3
#----------------- Space properties --------------------#
factor_size_cell = 1 # means : cell size = mean free path
# rectangle of size l1*l2
size_space = 3 # to be multiply by the mean free path
l1,l2 = size_space*mean_free_path,size_space*mean_free_path
# resolution along each previous dimension
res1, res2 = int(size_space/factor_size_cell), int(size_space/factor_size_cell)
nb_cells = res1*res2

if(debug): print("Dimension : {}x{} m".format(round(l1,2),round(l2,2)))

l3 = 0.1 # 10 cm => I use the deepness we shall use in practice.

#----------------- Grid creation ----------------------#
dtype = "LinkedList"
my_grid = Grid(l1,l2,[res1,res2], dtype = dtype) # LinkedList,DynamicArray

#-------------- Particles properties -------------------#
IODINE_MASS = get_mass_part(53, 53, 88)

charge = 0
mass = IODINE_MASS

# pos = MyVector(0,0,0), speed = MyVector(0,0,0), 
part_type = "I"
radius = effective_diameter/2.0
verbose = False
# status = -1 # irrelevant

#--------------- Particles creation -------------------#

# /!\ explanation of the calculus of Ne /!\
"""
Density of I : [I] = n = 10e20
size of the system : l1*l2*l3 (here)
Number of particles in the real system :  [I2]*l1*l2*l3 (if we inject them, we can go straight to here with the number of particles injected)
Remember that : mean free path is 0.01 m. (𝜆=1/(√2𝜋𝜎²𝑛)) (𝜎 : effective diameter of the particles)
At max speed (no electromagnetic field, no acceleration) : 3400 m/s (Note : Maxwellian distribution.)
That yields here : 𝑀𝑒𝑎𝑛𝐹𝑟𝑒𝑒𝑇𝑖𝑚𝑒 ≈ 3×10−6 𝑠 which yields Δ ≤ 3×10−7 𝑠. Which is what we'll choose.
In addition, we have a particles-per-cell-target of : Nc_mean = 50 particles
Which yields : Number_of_cells = res1*res2*Nc_mean particles.
"""

N_particles_real = int(real_particle_density*l1*l2*l3) # this is the REAL number of particles
N_particles_simu = int(mean_particles_number_per_cell*res1*res2)
Ne = int(N_particles_real/N_particles_simu)

if(debug): print("There is {} particles in the simulation. One particle accounts for {} real particles.".format("{:e}".format(N_particles_simu),"{:e}".format(Ne)))

def get_sigma(T, mass):
    # to be used for both sigma and maxwellian distribution
    # https://scicomp.stackexchange.com/questions/19969/how-do-i-generate-maxwell-boltzmann-variates-using-a-uniform-distribution-random
    # T in Kelvin
    kb = 1.38064e-23 # J K−1
    return np.sqrt(kb*T/mass) # constant mass for all particles

def get_T(v_mean, mass):
    kb = 1.38064e-23 # J K−1
    return mass*v_mean*v_mean/(3*kb)
sigma = get_sigma(T = get_T(300, mass), mass = mass) # T=454 K roughly for v_mean = 300 m/s
mu = 2*sigma*np.sqrt(2/np.pi) # only for maxwellian distribution

types = ['I']
numbers = [N_particles_simu]
speed_init_type = ['uniform'] # gaussian, uniform, uniform_norm, maxwellian : each axis speed will be initialized according to speed_init_type
speed_init_params = [[-300,300]]  # [[0,sigma]] # [[mu,sigma]], [[-300,300]], [[250,350]] # with if gaussian / maxwellian : mu = speed_init_params[0], sigma = speed_init_params[1], 
# if uniform (or uniform_norm) : min_speed = speed_init_params[0], max_speed = speed_init_params[1]
list_particles, vmean_exp = get_particles(types, numbers, speed_init_type, speed_init_params, effective_diameter, None, [0, 0], [l1,l2], verbose = False, debug = False)
print(get_T(300, mass))
print(mass)
#--------------- Rectangle creation -------------------#

rectangle_walls = [segment(Point(0,0),Point(0,l2)), segment(Point(0,0),Point(l1,0)), \
    segment(Point(l1,0),Point(l1,l2)), segment(Point(0,l2),Point(l1,l2))]

#---------- Adding particles to the grid -------------#
for particle in list_particles:
    my_grid.add(particle)

#---------- Creating collision handler ---------------#

    # update function
def f(Y,t,m,q):
    vx=Y[3]
    vy=Y[4]
    vz=Y[5]
    
    ax = 0 
    ay = 0
    az = 0
    return np.array([vx, vy, vz, ax, ay, az])

    # parameters
eta = 0
p = 0

    # DSMC param
DSMC_params = {
    'vr_max' : 2*mean_speed,
    'effective_diameter':  effective_diameter,
    'Ne' : Ne, # this is the number of real particles one simulated particle represents.
    'cell_volume' : l1*l2*l3/(res1*res2), # since we are in 2D I don't really know what to add here actually... For now, I add the 3rd dimension rough size, that is l3
    'mean_particles_number_per_cell': mean_particles_number_per_cell,
}
use_particles_collisions = False
use_DSMC = True
integration_scheme = euler_explicit # scipy_integrate_solve_ivp , rk4 , euler_explicit
collisionHandler = CollisionHandler(list_particles, rectangle_walls, f, eta, p, use_particles_collisions = use_particles_collisions, \
        use_DSMC = use_DSMC, grid = my_grid, DSMC_params = DSMC_params, integration_scheme=integration_scheme)

#------------------- Simulation -------------------#
if(save_test):
    params_dict = {
        'id_test' : str(id_test),
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
        'saving_period' : saving_period,
        'number_of_collisions' : 0,
        'mean_acceptance_rate' : 0,
        'mean_proba' : 0,
        'mean_vr_norm' : 0,
        'vr_max_init' : DSMC_params['vr_max'],
        'vr_max' : 0,
        'volume' : l1*l2*l3,
        'v_mean' : vmean_exp, 
    }
    data_analyser = DataSaver(list_particles, name_test = str(id_test), saving_directory = saving_directory)
    data_analyser.save_test_params(tests_summary_file_name, params_dict, use_saving_directory = False)
    # integration params
t = 0

if(debug): print("\nSTARTING SIMULATION...\n")
# simulation

elapsed_time = time()

if(save_test):
        data_analyser.save_everything_to_one_csv(erase = True)

if(not test_config):
    for k in tqdm(range(MAX_INTEGRATION_STEP)):
        #if(debug): print("\nStep {} over {}...\n".format(k+1, MAX_INTEGRATION_STEP))
        collisionHandler.step(dt, t, [])
        t+=dt
        if(save_test and (k%saving_period==0 or k == MAX_INTEGRATION_STEP-1)): # we are saving last frame
            # TODO : add params with "callback functions"
            # that we would do each time (or initialize them before in some way)
            # if it requires some initializing (parameters to set etc.)
            data_analyser.save_everything_to_one_csv()
            print(collisionHandler.save_collisions_matrix(name = "test_"+str(id_test)+"_collision_matrix", iteration = k))
    
    if(debug): print("\nElapsed  time for {} iterations with {} particules and with {} data structure : {}".format(MAX_INTEGRATION_STEP, N_particles_simu, dtype, round(time()-elapsed_time,3)))
    
    number_of_collisions = collisionHandler.get_collisions_count()
    mean_acceptance_rate = collisionHandler.get_acceptance_rate()
    mean_vr_norm = collisionHandler.get_mean_vr_norm()
    vr_max_final = collisionHandler.get_vr_norm()
    mean_proba = collisionHandler.get_mean_proba()
    # saving again to the csv.
    params_dict['mean_proba']=mean_proba
    params_dict['number_of_collisions']=number_of_collisions
    params_dict['mean_acceptance_rate']=mean_acceptance_rate
    params_dict['mean_vr_norm']=mean_vr_norm
    params_dict['vr_max']=vr_max_final

    data_analyser.update_saved_params(tests_summary_file_name, params_dict, use_saving_directory = False)

    if(debug): 
        print("\nNumber of collisions : {}".format("{:e}".format(number_of_collisions)))
        print("\nMean acceptance rate : {}".format(mean_acceptance_rate))

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
    collisionHandler.step(dt, t, [])
    if(debug) : print("\t[OK]")
# my modules imports
from modules.collisions_handler import CollisionHandler
from modules.particules import Particule
from modules.utils import segment
from modules.vector import MyVector
from modules.grid import Grid
from modules.utils import get_mass_part

# import for analysis
from modules.data_analysis import DataSaver, DataAnalyser
# imports
from dolfin import Point
from random import random
import numpy as np
from scipy.stats import maxwell

# comparaison
from time import time
# ----------------- saving test to file ? --------------- #
saving_directory = 'tests/' # be careful, it won't check if this the directory is already created ...
tests_summary_file_name = "tests_summary"
save_test = True
# ----------------- id test -------------------- #
id_test = 0
# --------------- Default analysis ? ---------------- #
perform_default_analysis = False
#----------------- debug parameters --------------------#
debug = True
#----------------- physics properties --------------------#

real_particle_density = 1e20 # for I
effective_diameter = 4e-10 # roughly the diameter of I
max_speed = 4e4 # 4000 m/s we say ~> Maxwellian distribution

mean_free_path = 1/(np.sqrt(2)*np.pi*effective_diameter*effective_diameter*real_particle_density)
mean_free_time = mean_free_path/max_speed
dt = 0.25 * mean_free_time

mean_particles_number_per_cell = 50

MAX_INTEGRATION_STEP = 10
#----------------- Space properties --------------------#
# resolution along each previous dimension
res1, res2 = 10, 10
nb_cells = res1*res2
# rectangle of size l1*l2
l1,l2 = res1*mean_free_path,res2*mean_free_path

if(debug): print("Dimension : {}x{} m".format(round(l1,2),round(l2,2)))

l3 = 0.1 # 10 cm ?

#----------------- Grid creation ----------------------#
dtype = "DynamicArray"
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
# TODO : add maxwellian distribution for speed

N_particles_real = int(real_particle_density*l1*l2*l3) # this is the REAL number of particles
N_particles_simu = int(mean_particles_number_per_cell*res1*res2)
Ne = int(N_particles_real/N_particles_simu)

if(debug): print("There is {} particles in the simulation. One particle accounts for {} real particles.".format("{:e}".format(N_particles_simu),"{:e}".format(Ne)))


init_type = "maxwellian" # uniform

    # Maxwellian distribution parameters
# parameters of the maxwellian distribution
# https://stackoverflow.com/questions/63300833/maxwellian-distribution-in-python-scipy
σ = 300
μ = 3000
a = σ * np.sqrt(np.pi/(3.0*np.pi - 8.0))
m = 2.0*a*np.sqrt(2.0/np.pi)
loc = μ - m

    # uniform distribution parameters
min_speed_uniform_distribution = 2500
max_speed_uniform_distribution = 3500

list_particles=[]

for k in range(N_particles_simu):
    # norm of the speed
    if(init_type=='maxwellian'):
        norm_speed = float(maxwell.rvs(loc, a))
    else :
        norm_speed = min_speed_uniform_distribution+random()*\
            (max_speed_uniform_distribution-min_speed_uniform_distribution)
    # direction of the speed
    theta = random()*2*np.pi
    cTheta = float(np.cos(theta))
    sTheta = float(np.sin(theta))
    # position
    x, y = random(), random()
    while(x==0.0 or x==1.0): # avoiding walls
        x=random()
    while(y==0.0 or y==1.0):
        y=random()  
    list_particles.append(Particule(charge = charge, radius = radius, 
        mass = mass, part_type = part_type, \
            speed=MyVector(norm_speed*cTheta,norm_speed*sTheta,0), \
                pos=MyVector(l1*random(),l2*random(),0), \
                    verbose = verbose))

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
    return [vx, vy, vz, ax, ay, az]

    # parameters
eta = 0
p = 0

    # DSMC params

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

DSMC_params = {
    'vr_max' : 2*max_speed,
    'effective_diameter':  effective_diameter,
    'Ne' : Ne, # this is the number of real particles one simulated particle represents.
    'cell_volume' : l1*l2*l3/(res1*res2) # since we are in 2D I don't really know what to add here actually... For now, I add the 3rd dimension rough size, that is l3
}

collisionHandler = CollisionHandler(list_particles, rectangle_walls, f, eta, p, use_particles_collisions = False, \
        use_DSMC = True, grid = my_grid, DSMC_params = DSMC_params)

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
        'loss_charge_proba' : p
    }
    data_analyser = DataSaver(list_particles, name_test = str(id_test), saving_directory = saving_directory)
    data_analyser.save_test_params(tests_summary_file_name, params_dict, use_saving_directory = False)
    # integration params
t = 0

if(debug): print("\nSTARTING SIMULATION...\n")
# simulation


elapsed_time = time()

for k in range(MAX_INTEGRATION_STEP):
    if(debug): print("\nStep {} over {}...\n".format(k+1, MAX_INTEGRATION_STEP))
    collisionHandler.step(dt, t, [])
    t+=dt
    if(save_test):
        data_analyser.save_everything_to_one_csv()

if(debug): print("\nElapsed  time for {} iterations with {} particules and with {} data structure : {}".format(MAX_INTEGRATION_STEP, N_particles_simu, dtype, round(time()-elapsed_time,3)))


if(debug): print("\nNumber of collisions : {}".format("{:e}".format(collisionHandler.get_collisions_count())))


if(perform_default_analysis):
    data_analyser = DataAnalyser(tests_summary_file_name)
    data_analyser.draw_speed_norm_distribution(id_test)
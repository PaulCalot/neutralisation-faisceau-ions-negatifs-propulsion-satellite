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
tests_summary_file_name = "tests_summary_v4"
save_test = True
saving_period = 10
# ----------------- id test -------------------- #
id_test = 14
# --------------- Default analysis ? ---------------- #
perform_default_analysis = False
#----------------- debug parameters --------------------#
debug = False
test_config = False 
#----------------- physics properties --------------------#

real_particle_density = 1e20 # for I
effective_diameter = 4e-10 # roughly the diameter of I
max_speed = 20000 # 20000 m/s we say ~> Maxwellian distribution
# TODO investigate max_speed importance

mean_free_path = 1/(np.sqrt(2)*np.pi*effective_diameter*effective_diameter*real_particle_density)
mean_free_time = mean_free_path/max_speed
dt = 0.25 * mean_free_time

mean_particles_number_per_cell = 50

MAX_INTEGRATION_STEP = 300
#----------------- Space properties --------------------#
# resolution along each previous dimension
res1, res2 = 10, 10
nb_cells = res1*res2
# rectangle of size l1*l2
l1,l2 = res1*mean_free_path,res2*mean_free_path

if(debug): print("Dimension : {}x{} m".format(round(l1,2),round(l2,2)))

l3 = 0.1 # 10 cm ?

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


init_type ="uniform" # "maxwellian" # uniform
speed_init_type ="2" # 2 

# type 1 : with theta and norm
# type 2 : with each vx, vy initialized and then normalized

    # Maxwellian distribution parameters
# parameters of the maxwellian distribution
# https://stackoverflow.com/questions/63300833/maxwellian-distribution-in-python-scipy
σ = 300
μ = 3000
a = σ * np.sqrt(np.pi/(3.0*np.pi - 8.0))
m = 2.0*a*np.sqrt(2.0/np.pi)
loc = μ - m

    # uniform distribution parameters
min_speed_uniform_distribution = 1800
max_speed_uniform_distribution = 2500

list_particles=[]
if(debug) : mean = 0
for k in range(N_particles_simu):
    # TODO : Distribution à faire sur vx / vy
    # norm of the speed
    my_speed = 0

    if(init_type=='1'):
        if(init_type=='maxwellian'):
            norm_speed = float(maxwell.rvs(loc, a))
        else:
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
            
        my_speed = MyVector(norm_speed*cTheta,norm_speed*sTheta,0)

        list_particles.append(Particule(charge = charge, radius = radius, 
            mass = mass, part_type = part_type, \
                speed=my_speed, \
                    pos=MyVector(l1*x,l2*y,0), \
                        verbose = verbose))

    elif(init_type=='2'):
        if(init_type=='maxwellian'): # TODO
            # we seek a temperature average of 3000 (for example)
            # we know that this everage is given by sqrt(8kT/(Pi m)) 
            # we only need T as we know the mass of other parts
            # a = sqrt(KT/m)
            const_k = 8.314
            get_T = lambda target : target*target*np.pi*mass/(8*const_k)
            get_a = lambda target : np.sqrt(get_T(target)*const_k/mass)
            target = μ
            a_ = get_a(target)
            vx = maxwell.rvs(0, a_)
            vy = maxwell.rvs(0, a_)
            vz = 0
            my_speed = target*MyVector(vx,vy,0).normalize()
        else:
            """
            norm_speed = min_speed_uniform_distribution+random()*\
                (max_speed_uniform_distribution-min_speed_uniform_distribution)
            vx = (1-2*random())
            vy = (1-2*random())
            my_speed = norm_speed*MyVector(vx,vy,0).normalize()
            """
            norm_speed = min_speed_uniform_distribution+random()*\
                (max_speed_uniform_distribution-min_speed_uniform_distribution)
            vx =  np.sign(1-2*random())*norm_speed
            norm_speed = min_speed_uniform_distribution+random()*\
                (max_speed_uniform_distribution-min_speed_uniform_distribution)
            vy = np.sign(1-2*random())*norm_speed
            my_speed = MyVector(vx,vy,0)
        x, y = random(), random()
        while(x==0.0 or x==1.0): # avoiding walls
            x=random()
        while(y==0.0 or y==1.0):
            y=random()  
        list_particles.append(Particule(charge = charge, radius = radius, 
            mass = mass, part_type = part_type, \
                speed=my_speed, \
                    pos=MyVector(l1*x,l2*y,0), \
                        verbose = verbose))    
    elif(init_type=='3'):
        if(init_type=='maxwellian'): # TODO
            # we seek a temperature average of 3000 (for example)
            # we know that this everage is given by sqrt(8kT/(Pi m)) 
            # we only need T as we know the mass of other parts
            # a = sqrt(KT/m)
            if(debug): print("Type {} for maxwellian distribution is not defined. Type 2 is used.".format(init_type))
            const_k = 8.314
            get_T = lambda target : target*target*np.pi*mass/(8*const_k)
            get_a = lambda target : np.sqrt(get_T(target)*const_k/mass)
            target = μ
            a_ = get_a(target)
            vx = maxwell.rvs(0, a_)
            vy = maxwell.rvs(0, a_)
            vz = 0
            my_speed = target*MyVector(vx,vy,0).normalize()
        else:
            min_speed_uniform_distribution = -4000
            max_speed_uniform_distribution = 4000
            vx = min_speed_uniform_distribution+random()*\
                (max_speed_uniform_distribution-min_speed_uniform_distribution)
            vy = min_speed_uniform_distribution+random()*\
                (max_speed_uniform_distribution-min_speed_uniform_distribution)
            my_speed = MyVector(vx,vy,0)

        x, y = random(), random()
        while(x==0.0 or x==1.0): # avoiding walls
            x=random()
        while(y==0.0 or y==1.0):
            y=random()  
        list_particles.append(Particule(charge = charge, radius = radius, 
            mass = mass, part_type = part_type, \
                speed=my_speed, \
                    pos=MyVector(l1*x,l2*y,0), \
                        verbose = verbose))    

    if(debug): 
        print(my_speed)
        mean += my_speed.norm()

if(debug): print("Mean speed init : {} m/s".format(round(mean/N_particles_simu,2)))
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
collisionHandler = CollisionHandler(list_particles, rectangle_walls, f, eta, p, use_particles_collisions = use_particles_collisions, \
        use_DSMC = use_DSMC, grid = my_grid, DSMC_params = DSMC_params, integration_scheme=integration_scheme)

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
        'init_type' : init_type,
        'speed_init_type' : speed_init_type
    }
    data_analyser = DataSaver(list_particles, name_test = str(id_test), saving_directory = saving_directory)
    data_analyser.save_test_params(tests_summary_file_name, params_dict, use_saving_directory = False)
    # integration params
t = 0

if(debug): print("\nSTARTING SIMULATION...\n")
# simulation


elapsed_time = time()

if(save_test):
        data_analyser.save_everything_to_one_csv()

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

    if(debug): print("\nElapsed  time for {} iterations with {} particules and with {} data structure : {}".format(MAX_INTEGRATION_STEP, N_particles_simu, dtype, round(time()-elapsed_time,3)))


    if(debug): print("\nNumber of collisions : {}".format("{:e}".format(collisionHandler.get_collisions_count())))


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
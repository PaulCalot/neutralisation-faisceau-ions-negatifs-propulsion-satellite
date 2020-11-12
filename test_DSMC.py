# my modules imports
from modules.collisions_handler import CollisionHandler
from modules.particules import Particule
from modules.utils import segment
from modules.vector import MyVector
from modules.grid import Grid
from modules.utils import get_mass_part

# imports
from dolfin import Point
from random import random
import numpy as np
from scipy.stats import maxwell

#----------------- debug parameters --------------------#
debug = True

#----------------- physics property --------------------#

real_particle_density = 1e20 # for I2
effective_diameter = 4e-10 # roughly the diameter of I2
max_speed = 4e4 # 4000 m/s we say ~> Maxwellian distribution

mean_free_path = 1/(np.sqrt(2)*np.pi*effective_diameter*effective_diameter*real_particle_density)
mean_free_time = mean_free_path/max_speed
dt = 0.25 * mean_free_time

mean_particles_number_per_cell = 50

MAX_INTEGRATION_STEP = 10
#----------------- Space properties --------------------#
# resolution along each previous dimension
res1, res2 = 10, 10

# rectangle of size l1*l2
l1,l2 = res1*mean_free_path,res2*mean_free_path

if(debug): print("Dimension : {}x{} m".format(round(l1,2),round(l2,2)))

l3 = 0.1 # 10 cm ?

#----------------- Grid creation ----------------------#

my_grid = Grid(l1,l2,[res1,res2])

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

if(debug): print("There is {} particles in the simulation. One particle account for {} real particles.".format("{:e}".format(N_particles_simu),"{:e}".format(Ne)))

# parameters of the maxwellian distribution
# https://stackoverflow.com/questions/63300833/maxwellian-distribution-in-python-scipy
σ = 300
μ = 3000

a = σ * np.sqrt(np.pi/(3.0*np.pi - 8.0))

m = 2.0*a*np.sqrt(2.0/np.pi)

loc = μ - m

list_particles=[]
for k in range(N_particles_simu):
    theta = random()*2*np.pi
    norm_speed = float(maxwell.rvs(loc, a))
    cTheta = float(np.cos(theta))
    sTheta = float(np.sin(theta))
    
    x = max(effective_diameter, min(random(),1.0-effective_diameter))
    y = max(effective_diameter, min(random(),1.0-effective_diameter))

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

# integration params
t = 0

if(debug): print("\nSTARTING SIMULATION...\n")
# simulation
for k in range(MAX_INTEGRATION_STEP):
    if(debug): print("\nStep {} over {}...\n".format(k+1, MAX_INTEGRATION_STEP))
    collisionHandler.step(dt, t, [])
    t+=dt
    #if(debug):
    #    for i, part in enumerate(list_particles):
    #       print("[{}] - pos : {} - speed : {}".format(i,(round(part.get_2D_pos().x,3),round(part.get_2D_pos().y,3)),\
    #           (round(part.get_2D_speed().x,3),round(part.get_2D_speed().y,3))))

if(debug): print("\nNumber of collisions : {}".format("{:e}".format(collisionHandler.get_collisions_count())))
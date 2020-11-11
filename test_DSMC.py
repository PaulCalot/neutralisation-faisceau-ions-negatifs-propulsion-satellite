# my modules imports
from modules.collisions_handler import CollisionHandler
from modules.particules import Particule
from modules.utils import segment
from modules.vector import MyVector
from modules.grid import Grid

# imports
from dolfin import Point
from random import random

#----------------- debug parameters --------------------#
debug = True

#----------------- Space properties --------------------#

# rectangle of size l1*l2
l1,l2 = 100.0,100.0

# resolution along each previous dimension
res1, res2 = 100, 100

#----------------- Grid creationg ----------------------#

my_grid = Grid(l1,l2,[res1,res2])

#-------------- Particles properties -------------------#

charge = 0
mass = 1
# pos = MyVector(0,0,0), speed = MyVector(0,0,0), 
part_type = "I"
radius = 1
verbose = False
# status = -1 # irrelevant

N_particles = 10

#--------------- Particles creation -------------------#
# TODO : add maxellian repartiton for speed
list_particles = [Particule(charge = charge, radius = radius, mass = mass, \
    part_type = part_type, speed=MyVector(l1*random(),l2*random(),0), \
        pos=MyVector(l1*(0.25+random()/2),l2*(0.25+random()/2),0), \
            verbose = verbose) for k in range(N_particles)]

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


collisionHandler = CollisionHandler(list_particles, rectangle_walls, f, eta, p)

#------------------- Simulation -------------------#

# integration params
dt = 0.01
t = 0
max_t = 20

# simulation
while t < max_t :
    collisionHandler.step(dt, t, [])
    t+=dt
    if(debug):
        for i, part in enumerate(list_particles):
            print("[{}] - pos : {} - speed : {}".format(i,(round(part.get_2D_pos().x,3),round(part.get_2D_pos().y,3)),\
                (round(part.get_2D_speed().x,3),round(part.get_2D_speed().y,3))))


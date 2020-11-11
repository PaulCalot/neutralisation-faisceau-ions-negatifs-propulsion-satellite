# TODO : this script does not work yet. I should fix the saving to mpeg stuff.

# my modules imports
from ..modules.collisions_handler import CollisionHandler
from .modules.particules import Particule
from .modules.utils import segment
from .modules.vector import MyVector

# imports
import matplotlib
matplotlib.use("Agg")
from dolfin import Point
from random import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# debugging
debug = True

# particules
parts_nb = 100
particules = [Particule(speed=MyVector(random(),random(),0),pos=MyVector(0.25+random()/2,0.25+random()/2,0)) for k in range(parts_nb)]

# walls
walls = [segment(Point(0,0),Point(0,1)), segment(Point(0,0),Point(1,0)), \
    segment(Point(1,0),Point(1,1)), segment(Point(0,1),Point(1,1))]

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

# collision handler
collisionHandler = CollisionHandler(particules, walls, f, eta, p)

# integration params
dt = 0.01
t = 0
max_t = 20

"""while t < max_t :
    collisionHandler.step(dt, t, [])
    t+=dt
    if(debug):
        for i,part in enumerate(particules):
            print("[{}] - pos : {} - speed : {}".format(i,(round(part.get_2D_pos().x,3),round(part.get_2D_pos().y,3)),\
                (round(part.get_2D_speed().x,3),round(part.get_2D_speed().y,3))))
"""

# For animation - works only in notebook
#%matplotlib notebook
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Set up formatting for the movie files - trying to save to mpeg
Writer = animation.FFMpegWriter()

writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
line, = ax.plot([], [], 'ro', ms=0.5)
ax.set_xlim(0,1)
ax.set_ylim(0,1)

# animation function
def animate(i):
    t = dt*i
    collisionHandler.step(dt, t, [])
    
    X_pos = []
    Y_pos = []
    for k in range(len(particules)):
        X_pos.append(particules[k].get_2D_pos().x)
        Y_pos.append(particules[k].get_2D_pos().y)
    line.set_data(X_pos, Y_pos)
    return line,



anim = animation.FuncAnimation(fig, animate, interval=max_t, blit=False)



anim.save('lines.mp4', writer=writer)
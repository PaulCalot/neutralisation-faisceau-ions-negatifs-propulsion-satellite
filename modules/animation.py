%matplotlib notebook
import matplotlib.pyplot as plt
import matplotlib.animation as animation

disks = Disks(15, 0.02)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
line, = ax.plot([], [], 'ro', ms=14)
ax.set_xlim(0,1)
ax.set_ylim(0,1)

def animate(i):
    disks.step(0.01)
    line.set_data(disks.positions[:,0], disks.positions[:,1])
    return line,
    
anim = animation.FuncAnimation(fig, animate, interval=10, blit=False)
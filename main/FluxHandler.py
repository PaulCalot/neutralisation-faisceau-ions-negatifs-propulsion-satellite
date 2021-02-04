from main import MyVector
from main import init_particles_flux

import matplotlib.pyplot as plt
import numpy as np

class FluxHandler(): # Maxwellian flux for now
    def __init__(self, list_walls, flux_in_wall, flux_out_wall, fluxes, \
        particles_types, speed_param1, speed_param2):
        self.particles_types=particles_types
        self.speed_params1 = speed_param1
        self.speed_params2 = speed_param2
        self.list_walls = list_walls 
        self.fluxes = fluxes
        self.flux_params = self.init_flux_params(flux_in_wall, flux_out_wall)

    def step(self, dt):
        nb_particles_to_inject = [int(flux * dt) for flux in self.fluxes]
        return init_particles_flux(self.flux_params['in_wall'], self.flux_params['in_direction'], nb_particles_to_inject, \
            self.particles_types, self.speed_param1, self.speed_param2)

    # ----------- init function -------------- #
    def select_wall_(self, criteria):
        # [x1, y1, x2, y2] = a wall
        wall = None
        direction = None # direction towards the system
        if(criteria == 'top'):
            wall = sorted(self.walls, key = lambda x: x[1]+x[3])[-1]
            direction = MyVector(wall[1]-wall[3],wall[2]-wall[0]).normalize()
            if(MyVector(0,1).inner(direction)>0):
                direction = -1.0*direction
        elif(criteria == 'bottom'):
            wall = sorted(self.walls, key = lambda x: x[1]+x[3])[0]
            direction = MyVector(wall[1]-wall[3],wall[2]-wall[0]).normalize()
            if(MyVector(0,1).inner(direction)<0):
                direction = -1.0*direction
        elif(criteria == 'left'):
            wall = sorted(self.walls, key = lambda x: x[0]+x[2])[0]
            direction = MyVector(wall[1]-wall[3],wall[2]-wall[0]).normalize()
            if(MyVector(1,0).inner(direction)<0):
                direction = -1.0*direction
        elif(criteria == 'right'):
            wall = sorted(self.walls, key = lambda x: x[0]+x[2])[-1]
            direction = MyVector(wall[1]-wall[3],wall[2]-wall[0]).normalize()
            if(MyVector(1,0).inner(direction)>0):
                direction = -1.0*direction
        return wall, direction

    def init_flux_params(self, in_flux, out_flux):
        in_direction = None
        out_direction = None
        if(in_flux == None):
            return None
        else:
            in_wall, in_direction = self.select_wall_(criteria = in_flux)
        if(out_flux == None):
            out_wall = None
        else:
            out_wall, out_direction = self.select_wall_(criteria = out_flux)
            out_direction = -1.0*out_direction
        return {
            'in_wall':in_wall,
            'in_direction': in_direction,
            'out_wall':out_wall,
            'out_direction':out_direction
        }

    def plot(self):
        if(self.flux != None):
            in_wall = self.flux['in_wall']
            if(in_wall != None):
                x, y = 0.5*(in_wall[0]+in_wall[2]), 0.5*(in_wall[1]+in_wall[3])
                size = np.sqrt((in_wall[2]-in_wall[0])**2+(in_wall[3]-in_wall[1])**2)
                in_direction = self.flux['in_direction']
                factor = size*0.2
                plt.arrow(x = x, y = y, dx = factor*in_direction.x, dy = factor*in_direction.y, color = 'r', width = 0.5*factor)
            out_wall = self.flux['out_wall']
            if(out_wall != None):
                x, y = 0.5*(out_wall[0]+out_wall[2]), 0.5*(out_wall[1]+out_wall[3])
                size = np.sqrt((out_wall[2]-out_wall[0])**2+(out_wall[3]-out_wall[1])**2)
                out_direction = self.flux['out_direction']
                factor =  size*0.2
                plt.arrow(x = x, y = y, dx = factor*out_direction.x, dy = factor*out_direction.y, color = 'r', width = 0.5*factor)
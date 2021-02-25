from main import MyVector
from main import init_particles_flux

import matplotlib.pyplot as plt
import numpy as np

class FluxHandler(): # Maxwellian flux for now
    debug=False
    def __init__(self):
        self.use_flux = False

    def step(self, dt):
        if(self.use_flux):
            nb_particles_to_inject_float = dt*self.fluxes
            remainders_, nb_particles_to_
            inject = np.modf(np.add(self.remainders, nb_particles_to_inject_float))
            nb_particles_to_inject = nb_particles_to_inject. astype(int)
            self.remainders = remainders_
            if(self.debug):
                
                print('Particles injected : {}, remainder : {}'.format(nb_particles_to_inject, self.remainders))
            return init_particles_flux(self.in_flux_wall, self.direction, nb_particles_to_inject, \
                self.particles_types, self.temperatures, self.drifts, dt)

    # ----------- init function -------------- #

    def init(self, list_walls, in_flux_wall, direction, fluxes, \
        particles_types, temperatures, drifts) -> None:
        self.particles_types=particles_types
        self.temperatures = temperatures
        self.drifts = drifts
        self.list_walls = list_walls 
        self.fluxes = np.array(fluxes)
        self.remainders = np.zeros(self.fluxes.shape)
        self.in_flux_wall = in_flux_wall
        self.direction = direction
        self.use_flux = True
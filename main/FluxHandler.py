from main import MyVector
from main import init_particles_flux

import matplotlib.pyplot as plt
import numpy as np

class FluxHandler(): # Maxwellian flux for now

    def __init__(self):
        self.use_flux = False

    def step(self, dt):
        if(self.use_flux):
            nb_particles_to_inject = [int(flux * dt) for flux in self.fluxes]
            return init_particles_flux(self.in_flux_wall, self.direction, nb_particles_to_inject, \
                self.particles_types, self.speed_param1, self.speed_param2)

    # ----------- init function -------------- #

    def init(self, list_walls, in_flux_wall, direction, fluxes, \
        particles_types, speed_param1, speed_param2) -> None:
        self.particles_types=particles_types
        self.speed_param1 = speed_param1
        self.speed_param2 = speed_param2
        self.list_walls = list_walls 
        self.fluxes = fluxes
        self.in_flux_wall = in_flux_wall
        self.direction = direction
        self.use_flux = True
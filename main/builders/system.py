from main import FluxHandler
from main import available_particles
from main import get_maxwellian_mean_speed_from_temperature


import numpy as np
from abc import ABC, abstractmethod
from main import DSMC

# plotting 
import matplotlib.pyplot as plt

from main import MyVector

class system(ABC):
    debug = False

    def __init__(self, options, min_mean_free_path) -> None:
        self.options = options
        self.min_mean_free_path = min_mean_free_path
        self.size = [0,0,0]
        self.resolution = [0,0]
        self.offset = [0,0]
        # initializing - will call subclasses function if called
        # first walls then grid
        self.walls = self.make_system()
        #self.update_walls_with_offset()
        self.grid = self.make_grid()
        # other inits
        self.volume = self.init_system_volume()
        self.Ne = self.init_Ne()
        
        # init flux if there is
        # TODO : init the flux how it should be!
        # flux_params : 
            # 'in_wall':in_wall,
            # 'in_direction': in_direction,
            # 'out_wall':out_wall,
            # 'out_direction':out_direction
        self.temperatures = self.options['temperatures']
        self.flux_params = self.init_flux_params(self.options['flux_in_wall'], self.options['flux_out_wall'])
        self.flux_handler = FluxHandler()
        if(self.flux_params != None):
            l1, l2 = len(self.options['particles_types']), len(self.temperatures)
            self.temperatures = [self.temperatures[0]]*l1 if(l2!=l1) else self.temperatures #
            self.fluxes = self.init_fluxes() # consist of initializing gamma which depends upon particles density etc.
            self.flux_handler.init(self.walls, self.flux_params['in_wall'], self.flux_params['in_direction'], \
                self.fluxes, self.options['particles_types'], self.options['temperatures'], self.options['drifts'])
            if(self.debug):
                print("Fluxes : {}".format(self.fluxes))
    
        self.factor = .5*(self.size[0]/self.resolution[0]+self.size[1]/self.resolution[1])

        # added for dictionnary purposes 
        self.mean_speed = None

        if(self.debug):
            print('Volume : {} vs {} m-3'.format(self.volume, self.size[0]*self.size[1]*self.size[2]))
            print('Ne : {} '.format(self.Ne))
    # the abstracts methods mean they need to be redefined by the subclass
    @abstractmethod
    def make_grid(self):
        pass

    @abstractmethod
    def make_system(self):
        pass

    def init_system_volume(self):
        return self.size[2]*self.grid.get_surface()

    def update_walls_with_offset(self):
        ox, oy = self.offset
        walls = []
        for wall in self.walls:
            x1,y1,x2,y2 = wall
            wall = [x1+ox,y1+oy, x2+ox, y2+oy]
            walls.append(wall)
        self.walls = walls

    def init_Ne(self):
        if(self.debug):
            print('Number of cells: {} vs {}'.format(self.grid.get_number_of_cells(), self.resolution[0]*self.resolution[1]))
        N_particles_real = [int(density*self.volume) for density in self.options['particles_densities']] # this is the REAL number of particles
        N_particles_simu = [int(nb*self.grid.get_number_of_cells()) for nb in self.options['particles_mean_number_per_cell']]
        return [i/j for i, j in zip(N_particles_real,N_particles_simu)]

    def plot(self):
        self.grid.plot()
        # adding now the walls
        fig = plt.plot(figsize = (10,30))
        for wall in self.walls :
            self.plot_wall(wall)
        self.plot_fluxes()
        plt.axis('scaled')
        #plt.savefig('thruster_with_arrows.png', dpi = 1200)

    def add(self, list_particles): # should not be use
        self.dsmc.add_particles(list_particles)
    
    def get_size(self):
        return self.size
        
    def get_resolutions(self):
        return self.resolution
    
    def get_walls(self):
        return self.walls
    
    def get_grid(self):
        return self.grid
    
    def get_resolution(self):
        return self.resolution

    def get_offset(self):
        return self.offset

    def plot_wall(self,wall):
        # wall = [x1,y1,x2,y2]
        x1,y1,x2,y2 = wall
        plt.plot([x1,x2], [y1,y2], color = 'k')

    def get_number_of_particles(self):
        nb = self.grid.get_number_of_particles() if self.flux_params == None else self.dsmc.get_number_of_particles()
        return nb
    
    def get_number_of_cells(self):
        return self.grid.get_number_of_cells()
        
    def get_list_particles(self):
        return self.dsmc.get_list_particles()

    def get_mean_speed(self):
        return self.mean_speed
    
    def get_integration_scheme(self):
        return self.dsmc.get_scheme()

    def get_volume(self):
        return self.volume
    
    def get_Ne(self):
        return self.Ne
        
    # ---------- step function ------------ #
    
    def step(self, dt, t):
        self.dsmc.step(dt, t, self.f_args)
        new_particles = self.flux_handler.step(dt)
        if(self.debug):print('Injecting {} new particles.'.format(len(new_particles)))
        self.dsmc.add_particles(new_particles)
        
    def simulate(self, dt, max_number_of_steps):
        return
    # ----------- init DSMC --------------- #
    def init_dsmc(self, mean_speed, integration_scheme, f, f_args) -> None:
        # TODO : make it better
        self.f_args =f_args
        effective_diameters = [available_particles[type_]['effective diameter'] for type_ in self.options['particles_types']]
        lx, ly, lz = self.size 
        res_x, res_y = self.resolution
        if(mean_speed==None):
            masses = [available_particles[type_]['mass'] for type_ in self.options['particles_types']]
            v_mean_list = np.array([get_maxwellian_mean_speed_from_temperature(self.temperatures[k], masses[k]) for k in range(len(masses))])
            mean_speed = np.mean(v_mean_list)
        self.mean_speed = mean_speed
        dsmc_params = {
            'vr_max' : 2*mean_speed,
            'effective_diameter':  effective_diameters[0], # for now
            'Ne' : self.Ne[0], # this is the number of real particles one simulated particle represents.
            'cell_volume' : lx*ly*lz/(res_x*res_y), # since we are in 2D I don't really know what to add here actually... For now, I add the 3rd dimension rough size, that is l3
            'mean_particles_number_per_cell': sum(self.options['particles_mean_number_per_cell']) # we don't account for the type
        }
        out_wall = self.flux_params['out_wall'] if self.flux_params != None else None
        self.dsmc = DSMC(self, dsmc_params, integration_scheme, f, out_wall = out_wall)

    # ----------- out / in flux params ------------- #
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
            in_wall = [in_wall[0],in_wall[1],in_wall[2],in_wall[3]]
        if(out_flux == None):
            out_wall = None
        else:
            out_wall, out_direction = self.select_wall_(criteria = out_flux)
            out_wall = [out_wall[0],out_wall[1],out_wall[2],out_wall[3]]
            out_direction = -1.0*out_direction
        return {
            'in_wall':in_wall,
            'in_direction': in_direction,
            'out_wall':out_wall,
            'out_direction':out_direction
        }
    
    def init_fluxes(self):
        masses = [available_particles[type_]['mass'] for type_ in self.options['particles_types']]
        v_mean_list = [get_maxwellian_mean_speed_from_temperature(self.temperatures[k], masses[k]) for k in range(len(masses))]
        in_wall = self.flux_params['in_wall']
        lengh_wall = np.sqrt((in_wall[2]-in_wall[0])**2+(in_wall[3]-in_wall[1])**2)
        if(self.debug):
            print('Lenght wall : {}'.format(lengh_wall))
        surface_in = lengh_wall*self.size[2] # lengh_wall*lz
        densities = self.options['particles_densities']

        # function
        get_flux = lambda n, v_mean, drift : 0.25*n*v_mean*surface_in+n*drift*surface_in # this function should be modified to accomodate the true flux

        drifts = self.options['drifts']
        fluxes = [get_flux(n, v_mean, drift) for n, v_mean, drift in zip(densities, v_mean_list, drifts)]

        fluxes = [flux/N for flux, N in zip(fluxes, self.Ne)]
        if(self.debug):print('My fluxes {}'.format(fluxes))
        return fluxes # particles/s
    
    def plot_fluxes(self):
        if(self.flux_params != None):
            in_wall = self.flux_params['in_wall']
            if(in_wall != None):
                x, y = 0.5*(in_wall[0]+in_wall[2]), 0.5*(in_wall[1]+in_wall[3])
                size = np.sqrt((in_wall[2]-in_wall[0])**2+(in_wall[3]-in_wall[1])**2)
                in_direction = self.flux_params['in_direction']
                factor = size*0.2
                plt.arrow(x = x, y = y, dx = factor*in_direction.x, dy = factor*in_direction.y, color = 'r', width = 0.5*factor)
            out_wall = self.flux_params['out_wall']
            if(out_wall != None):
                x, y = 0.5*(out_wall[0]+out_wall[2]), 0.5*(out_wall[1]+out_wall[3])
                size = np.sqrt((out_wall[2]-out_wall[0])**2+(out_wall[3]-out_wall[1])**2)
                out_direction = self.flux_params['out_direction']
                factor =  size*0.2
                plt.arrow(x = x, y = y, dx = factor*out_direction.x, dy = factor*out_direction.y, color = 'r', width = 0.5*factor)
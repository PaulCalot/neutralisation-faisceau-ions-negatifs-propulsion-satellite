import numpy as np
from main.data_structures.grid import Grid

class system:
    def __init__(self, options, min_mean_free_path) -> None:
        self.options = options
        self.min_mean_free_path = min_mean_free_path
        self.factor, self.size = 0, [0,0]
        self.grid = self.make_grid()
        self.walls = self.make_walls() # TBD by the child class 
        self.size = [0,0,0]

    def make_grid(self):
        self.factor = self.options['factor']
        self.resolution = self.options['size']
        return Grid(self.resolution[0]*self.factor,self.resolution[1]*\
            self.factor,[self.resolution[0],self.resolution[1]], dtype = 'LinkedList')
    
    def make_walls(self):
        return

    def plot(self):
        self.grid.plot()

    def add(self, particle):
        self.grid.add(particle)

    def get_size(self):
        return self.size
        
    def get_space_resolutions(self):
        return self.resolution
    
    def get_walls(self):

        return self.walls
    
    def get_grid(self):
        return self.grid
    


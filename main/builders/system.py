import numpy as np
from abc import ABC, abstractmethod

# plotting 
import matplotlib.pyplot as plt

class system(ABC):
    def __init__(self, options, min_mean_free_path) -> None:
        self.options = options
        self.min_mean_free_path = min_mean_free_path
        self.factor, self.size = 0, [0,0]
        self.size = [0,0,0]
        self.resolution = [0,0]
        self.offset = None

        # initializing - will call subclasses function if called
        self.walls = self.make_system()
        self.grid = self.make_grid()
        
    # the abstracts methods mean they need to be redefined by the subclass
    @abstractmethod
    def make_grid(self):
        pass

    @abstractmethod
    def make_system(self):
        pass

    def plot(self):
        self.grid.plot()
        # adding now the walls 
        for wall in self.walls :
            self.plot_wall(wall)
        
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
    
    def get_resolution(self):
        return self.resolution

    def get_offset(self):
        return self.offset

    def plot_wall(self,wall):
        # wall = [x1,y1,x2,y2]
        x1,y1,x2,y2 = wall
        plt.plot([x1,x2], [y1,y2], color = 'k')

from .system import system
from main import segment
from dolfin import Point
from main.data_structures.grid import Grid
import warnings

class square(system):
    def __init__(self, options, min_mean_free_path = None) -> None:
        super().__init__(options, min_mean_free_path)
    
    def make_system(self): # called first
        self.size = self.options['size']
        self.resolution = self.options['resolution']
        lx, ly = self.size[0], self.size[1]
        cell_size_max = max(lx/self.resolution[0], ly/self.resolution[1])
        if(self.min_mean_free_path != None and cell_size_max>self.min_mean_free_path):
            message = "Max cell size is {:e} m and higher than the min mean free path which is {:e}.".format(cell_size_max,self.min_mean_free_path)
            warnings.warn(message)
        self.size = [lx,ly,self.options['lz']]
        walls = [segment(Point(0,0),Point(0,ly)), segment(Point(0,0),Point(lx,0)), \
            segment(Point(lx,0),Point(lx,ly)), segment(Point(0,ly),Point(lx,ly))]
        return walls

    def make_grid(self): # called second
        return Grid(self.size[0],self.size[1],[self.resolution[0],self.resolution[1]], dtype = 'LinkedList')

    
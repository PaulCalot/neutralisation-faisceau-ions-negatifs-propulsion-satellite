from .system import system
from main import segment
from dolfin import Point
from main.data_structures.grid import Grid

class square(system):
    def __init__(self, options, min_mean_free_path) -> None:
        super().__init__(options, min_mean_free_path)
    
    def make_grid(self):
        return Grid(self.size[0],self.size[1],[self.resolution[0],self.resolution[1]], dtype = 'LinkedList')

    def make_system(self):
        self.factor = self.options['factor']
        self.resolution = self.options['size']
        lx, ly = self.min_mean_free_path*self.resolution[0]*self.factor,\
            self.min_mean_free_path*self.resolution[1]*self.factor
        self.size = [lx,ly,self.options['lz']]
        walls = [segment(Point(0,0),Point(0,ly)), segment(Point(0,0),Point(lx,0)), \
            segment(Point(lx,0),Point(lx,ly)), segment(Point(0,ly),Point(lx,ly))]
        return walls

    
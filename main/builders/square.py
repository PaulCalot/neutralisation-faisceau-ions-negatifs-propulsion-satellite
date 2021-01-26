from .system import system
from main import segment
from dolfin import Point

class square(system):
    def __init__(self, options, min_mean_free_path) -> None:
        super().__init__(options, min_mean_free_path)
        self.walls = self.make_walls()

    def make_walls(self):
        lx, ly = self.min_mean_free_path*self.resolution[0]*self.factor,\
            self.min_mean_free_path*self.resolution[1]*self.factor
        self.size = [lx,ly,self.options['lz']]
        walls = [segment(Point(0,0),Point(0,ly)), segment(Point(0,0),Point(lx,0)), \
            segment(Point(lx,0),Point(lx,ly)), segment(Point(0,ly),Point(lx,ly))]
        return walls

    
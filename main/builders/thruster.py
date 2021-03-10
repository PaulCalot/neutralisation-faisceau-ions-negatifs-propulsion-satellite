from .system import system
from main import segment, get_mesh, get_VandE
from dolfin import Point
from main.data_structures.grid import Grid
from main import init_particles_in_system
class thruster(system):
    debug = False
    def __init__(self, options, min_mean_free_path=None) -> None:
        super().__init__(options, min_mean_free_path)
        # mesh may not be needed right now but may be later.
        #print("Offset vs grid offset : {} vs {}".format(self.offset, self.grid.offsets))
        
    def make_system(self):
        mesh, segments_list, zone = get_mesh(self.options['geometry'])
        
        phi, E = get_VandE(mesh, self.options['geometry'], self.options['phi'])
        
        # rectangle of size l1*l2
        segment_X_list = []
        segment_Y_list = []
        for segment in segments_list:
            segment_X_list.append(segment[0])
            segment_X_list.append(segment[2])
            segment_Y_list.append(segment[1])
            segment_Y_list.append(segment[3])
        max_x, min_x = max(segment_X_list), min(segment_X_list)
        max_y, min_y = max(segment_Y_list), min(segment_Y_list)

        self.offset = [-min_x, -min_y]
        
        lx,ly = max_x - min_x, max_y - min_y
        lz = self.options['general']['lz']
        self.size = [lx,ly,lz] # offset not taken into account here. It's purely the size of the system
        # and thus it's purely positive...

        if(self.debug): print("Dimension : {}x{}x{} m".format(round(lx,2),round(ly,2),lz))        

        self.mesh = mesh
        self.zone = zone
        self.E = E
        self.phi = phi

        return segments_list

    def make_grid(self):
        lx,ly = self.size[0], self.size[1]
        self.resolution = self.options['general']['resolution']
        res_x, res_y = self.resolution[0], self.resolution[1]
        offset_x, offset_y = self.offset[0], self.offset[1]
        grid = Grid(lx,ly,[res_x,res_y], offsets = [offset_x,offset_y], dtype = 'LinkedList') # LinkedList , DynamicArray
        self.set_sparsed_space(grid)
        return grid

    def get_zone(self):
        return self.zone
    
    def set_sparsed_space(self, grid):
        particles_list = list(init_particles_in_system(['I'], [1e20], [100], ['gaussian'], [200],\
     [200], self.size, self.resolution, zone = self.zone, offsets = self.offset, verbose = False, debug = self.debug,)[0])
        grid.fill_sparsed_space_from_initial_particles_position(particles_list)

    def plot(self):
        super().plot()
        import matplotlib.pyplot as plt
        fig = plt.gcf()
        fig.set_size_inches(4, 20)
        #plt.rcParams["figure.figsize"] = (4,20)       
        plt.savefig('thruster_grid.png')
        #from fenics import plot
        #plot(self.phi)
        #from fenics import sqrt, dot
        #Ex, Ey = self.E.split(deepcopy=True)
        #NE=sqrt(dot(self.E,self.E))

    def get_zone(self):
        return self.zone
    
    def get_E(self):
        return self.E

    
from main import DynamicArray, LinkedList, MyVector

import numpy as np
import matplotlib.pyplot as plt

# 2 grid to start with
class Grid(object):
    # TODO : make it 3D (or rather nD)
    debug = False

    def __init__(self, lx, ly, resolutions, offsets = None, dtype = "LinkedList"):
        self.lx = lx
        self.ly = ly
        self.res = [resolutions, resolutions] if len(resolutions) == 1 else resolutions # can be res if it's the same for all directions, or [l_res, h_res]
        self.data_structure_class = DynamicArray if (dtype == "DynamicArray") else LinkedList
        self.use_offsets = False if offsets == None else True
        if(self.use_offsets):
            self.offsets = offsets if len(offsets)==2 else [offsets,offsets]
            self.offsets = MyVector(self.offsets[0], self.offsets[1], 0)
        else:
            self.offsets = MyVector(0,0,0)
        self.dtype = dtype
        self.grid = np.empty((self.res[0], self.res[1]), dtype = self.data_structure_class)          
        self.sparsed_space_idx = None
        self.number_of_particles = 0

    def add(self, particle):
        pos = particle.get_pos()+self.offsets
        if(self.debug) : print("Position with offsets correcting : {}".format(pos))
        self.add_(particle, self.get_pos_in_grid(pos))
        self.number_of_particles += 1

    def add_(self, particle, pos_in_grid):
        if(self.debug) : print("Resulting position in grid {}".format(pos_in_grid))
        if(self.debug):
            print("Adding ... " + particle.to_string() + " in grid position {}.".format(pos_in_grid), end=" " )
        
        if(self.grid[pos_in_grid[0],pos_in_grid[1]] == None): # checking if we already created this list or not
            self.grid[pos_in_grid[0],pos_in_grid[1]] = self.data_structure_class()
            if(self.debug):
                print("Created a new data structure of type {}.".format(self.dtype), end = " ")
        else:
            if(self.debug):
                print("Added the particule to pre-existing data structure (type {})".format(self.dtype), end = ' ')
        self.grid[pos_in_grid[0],pos_in_grid[1]].insert(particle)
        
        if(self.debug):
            print("     [OK]")

    def add_and_verify(self, particle):
        pos = particle.get_pos()+self.offsets
        if(self.debug) : print("Position with offsets correcting : {}".format(pos))
        pos_in_grid = self.get_pos_in_grid(pos)
        if(pos_in_grid == None):
            # Not in grid
            return False
        self.add_(particle, pos_in_grid)
        self.number_of_particles += 1

    def remove(self, particule):
        pos = particule.get_pos()+self.offsets
        self.remove_(particule, self.get_pos_in_grid(pos))

    def remove_(self, particule, pos_in_grid):
        if(self.debug): print("\nRemoving {} in position {}.".format(particule.to_string(), pos_in_grid))
        self.grid[pos_in_grid[0],pos_in_grid[1]].delete(particule)

    def remove_with_old_pos(self,particule,old_postion):
        pos = old_postion+self.offsets
        self.remove_(particule, self.get_pos_in_grid(pos))

    def update(self, particule, old_position):
        #if(self.debug):
        #    print("Updating position ... " + particule.to_string(), end= " ")

        pos = particule.get_pos()+self.offsets
        pos_in_grid = self.get_pos_in_grid(pos)
        
        if pos_in_grid == None : return False

        # useful for the system thruster amongst other...
        if(self.sparsed_space_idx != None and pos_in_grid not in self.sparsed_space_idx):
            return False

        old_pos_in_grid = self.get_pos_in_grid(old_position+self.offsets)
        # in theory old_pos_in_grid can not be out of the grid.

        if(self.debug): print("Current position in grid {} - Old position {}".format(pos_in_grid, old_pos_in_grid))

        if(pos_in_grid != old_pos_in_grid):
            if(self.debug):
                print("Positions in grid were different ...", end = "   [OK]")
            try :
                # may be there is some weird stuff going on and i should differentiate add_ and remove_
                self.add_(particule, pos_in_grid)
            except IndexError: # (ValueError,):
                print("New position {} not in grid.".format(pos_in_grid))
                return False
            self.remove_(particule, old_pos_in_grid)
        else :
            if(self.debug):
                print("Same positions.", end = "   [OK]")

        if(self.debug):
            print("     [OK]")
        return True

    def get_pos_in_grid(self, position):
        #return [int(position.x*self.res[0]/self.lx), int(position.y*self.res[1]/self.ly)]
        pos_x = position.x*self.res[0]/self.lx
        # First problem : we accepted -1 because list[-1] returns the last element which works. So I add to add the final if.
        # Second problem : there were a bug because int(-0.5) = 0 ... which works despite being out of the grid. I had to add the double if/else.
        if(pos_x < 0):
            pos_x = -1
        else :
            pos_x = int(pos_x)
        
        pos_y = position.y*self.res[1]/self.ly

        if(pos_y < 0):
            pos_y = -1
        else :
            pos_y = int(pos_y)

        
        if(self.debug) : print("Real position : [{},{}], position in grid : [{},{}].".format(position.x, position.y, pos_x, pos_y))
        #return [min(max(0,pos_x),self.res[0]-1), min(max(0,pos_y),self.res[1]-1)]
        if(pos_x >= self.res[0] or pos_x < 0 or pos_y >= self.res[1] or pos_y < 0):
            # not in the grid
            return None
        return [pos_x, pos_y]
    
    def get_closest_particules(self, particule, return_list = True):
        pos = particule.get_pos()+self.offsets
        pos_in_grid = self.get_pos_in_grid(pos)
        data_structure = self.grid[pos_in_grid[0],pos_in_grid[1]]
        if(return_list):
            return data_structure.to_list()
        else :
            return data_structure

    def list_to_string(self, position):
        pos_in_grid = self.get_pos_in_grid(position)
        self.list_to_string_(pos_in_grid)

    def list_to_string_(self, pos_in_grid):
        pos_list = self.grid[pos_in_grid[0],pos_in_grid[1]]
        if(pos_list != None):
            pos_list = pos_list.to_list()
            print('\nPosition in grid : {}'.format(pos_in_grid))
            for k in range(len(pos_list)):
                print(pos_list[k].to_string())
            
    def to_string(self):
        print("\nParticules in the grid : ")
        for i in range(self.res[0]):
            for j in range(self.res[1]):
                pos_in_grid = [i,j]
                self.list_to_string_(pos_in_grid)
    
    def get_all_particles(self):
        particles=[]
        for i in range(self.res[0]):
            for j in range(self.res[1]):
                data_structure = self.grid[i,j]
                if(data_structure!=None):
                    for part in self.grid[i,j].to_list():
                        particles.append(part)
        return particles

    def plot_cell_(self, i, j):
        dx, dy = self.lx/self.res[0],self.ly/self.res[1]
        x, y = self.lx*i/self.res[0]-self.offsets.x, self.ly*j/self.res[1]-self.offsets.y
        plt.plot([x, x+dx, x+dx, x, x], [y, y, y+dy, y+dy, y], color='g', alpha = 0.5)

    def plot_grid(self):
        # the goal here is to plot the grid 
        for i in range(self.res[0]):
            for j in range(self.res[1]):
                if(self.sparsed_space_idx == None or [i,j] in self.sparsed_space_idx):
                    self.plot_cell_(i,j)

    def plot(self):
        particles = self.get_all_particles()
        positions = [part.get_pos() for part in particles]
        speed = [part.get_speed() for part in particles]
        norm = [v.norm() for v in speed]
        x_pos, y_pos = [pos.x for pos in positions], [pos.y for pos in positions]

        plt.scatter(x_pos, y_pos, s=0.3, c = norm, cmap='seismic') # problem if offset
        self.plot_grid()

    # ------------ Getter and setter ------------- #

    def get_grid(self):
        return self.grid
    
    def get_number_of_particles(self):
        return self.number_of_particles

    def get_number_of_cells(self):
        if(self.sparsed_space_idx==None):
            return self.res[0]*self.res[1]
        else :
            return len(self.sparsed_space_idx)    
    
    def get_surface(self):
        # size of cell
        dx, dy = self.lx/self.res[0], self.ly/self.res[1]
        return dx*dy*self.get_number_of_cells()

    # ------------ Sparsed Space ---------------- #
    def fill_sparsed_space_from_initial_particles_position(self, list_particles = None):
        self.sparsed_space_idx = []
        list_particles = list_particles if list_particles != None else self.get_all_particles()
        for part in list_particles:
            pos = part.get_pos()+self.offsets
            pos_x = int(pos.x*self.res[0]/self.lx)
            pos_y = int(pos.y*self.res[1]/self.ly)
            self.sparsed_space_idx.append([pos_x,pos_y])

        # not optimized by not used  during the simulation so its okay
        temp = []
        for i in self.sparsed_space_idx:
            if i not in temp:
                temp.append(i)
        self.sparsed_space_idx = temp # remove same elements

        
       # print(self.sparsed_space_idx)

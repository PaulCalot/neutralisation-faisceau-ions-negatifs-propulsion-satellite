
import numpy as np

# local import
from .particules import Particule
from .dynamic_arrays import DynamicArray
from .linkedList import LinkedList 

# 2 grid to start with
class Grid(object):
    # TODO : make it 3D (or rather nD)
    debug = False

    def __init__(self, lx, ly, resolutions, dtype = "LinkedList"):
        self.lx = lx
        self.ly = ly
        self.res = [resolutions, resolutions] if len(resolutions) == 1 else resolutions # can be res if it's the same for all directions, or [l_res, h_res]
        self.data_structure_class = DynamicArray if (dtype == "DynamicArray") else LinkedList
        self.dtype = dtype
        self.grid = np.empty((self.res[0], self.res[1]), dtype = self.data_structure_class)

    def add(self, particule):
        pos = particule.get_pos()
        self.add_(particule, self.get_pos_in_grid(pos))

    def add_(self, particule, pos_in_grid):
    
        if(self.debug):
            print("Adding ... " + particule.to_string() + " in grid position {}.".format(pos_in_grid), end=" " )
        
        if(self.grid[pos_in_grid[0],pos_in_grid[1]] == None): # checking if we already created this list or not
            self.grid[pos_in_grid[0],pos_in_grid[1]] = self.data_structure_class()
            if(self.debug):
                print("Created a new data structure of type {}.".format(self.dtype), end = " ")
        else:
            if(self.debug):
                print("Added the particule to pre-existing data structure (type {})".format(self.dtype), end = ' ')
        self.grid[pos_in_grid[0],pos_in_grid[1]].insert(particule)
        
        if(self.debug):
            print("     [OK]")

    def remove(self, particule):
        pos = particule.get_pos()
        self.remove_(particule, self.get_pos_in_grid(pos))

    def remove_(self, particule, pos_in_grid):
        if(self.debug): print("\nRemoving {} in position {}.".format(particule.to_string(), pos_in_grid))
        self.grid[pos_in_grid[0],pos_in_grid[1]].delete(particule)

    def update(self, particule, old_position):
        if(self.debug):
            print("Updating position ... " + particule.to_string(), end= " ")

        pos = particule.get_pos()
        pos_in_grid = self.get_pos_in_grid(pos)
        old_pos_in_grid = self.get_pos_in_grid(old_position)

        if(pos_in_grid != old_pos_in_grid):
            if(self.debug):
                print("Positions in grid were different ...", end = "   [OK]")
            self.add_(particule, pos_in_grid)
            self.remove_(particule, old_pos_in_grid)
        if(self.debug):
            print("     [OK]")

    def get_pos_in_grid(self, position):
        #return [int(position.x*self.res[0]/self.lx), int(position.y*self.res[1]/self.ly)]
        pos_x = int(position.x*self.res[0]/self.lx)
        pos_y = int(position.y*self.res[1]/self.ly)
        return [min(max(0,pos_x),self.res[0]-1), min(max(0,pos_y),self.res[1]-1)]
        #return [pos_x, pos_y]

    # TODO : change theses functions to take into account both data structure
    # for now we don't use it so I ignore it.
    
    def get_closest_particules(self, particule, return_list = True):
        pos = particule.get_pos()
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
    
    # ------------ Getter and setter ------------- #

    def get_grid(self):
        return self.grid
    
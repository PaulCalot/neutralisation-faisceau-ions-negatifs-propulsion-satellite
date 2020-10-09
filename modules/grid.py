import numpy as np

# local import
from .particules import Particule

# 2 grid to start with
class Grid(object):
    # TODO : make it 3D (or rather nD)
    debug = False

    def __init__(self, lx, ly, resolutions):
        self.lx = lx
        self.ly = ly
        self.res = [resolutions, resolutions] if len(resolutions) == 1 else resolutions # can be res if it's the same for all directions, or [l_res, h_res]
        self.grid = np.empty((self.res[0], self.res[1]), dtype = LinkedList)

    def add(self, particule):
        pos = particule.get_pos()
        self.add_(particule, self.get_pos_in_grid(pos))

    def add_(self, particule, pos_in_grid):
    
        if(self.debug):
            print("Adding ... " + particule.to_string() + " in grid position {}.".format(pos_in_grid), end=" " )
        if(self.grid[pos_in_grid[0],pos_in_grid[1]] == None): # checking if we already created this list or not
            self.grid[pos_in_grid[0],pos_in_grid[1]] = LinkedList()
            if(self.debug):
                print("Created a new linkedlist.", end = " ")
        else:
            if(self.debug):
                print("Added the particule to pre-existing linkelist", end = ' ')
        self.grid[pos_in_grid[0],pos_in_grid[1]].insert(particule)
        if(self.debug):
            print("     [OK]")

    def remove(self, particule):
        pos = particule.get_pos()
        self.remove_(particule, self.get_pos_in_grid(pos))

    def remove_(self, particule, pos_in_grid):
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
        return [int(position.x*self.res[0]/self.lx), int(position.y*self.res[1]/self.ly)]

    def get_closest_particules(self, particule, return_list = True):
        pos = particule.get_pos()
        pos_in_grid = self.get_pos_in_grid(pos)
        linkedlist = self.grid[pos_in_grid[0],pos_in_grid[1]]
        if(return_list):
            list_ = []
            current = linkedlist.get_head()
            while(current):
                list_.append(current.get_data())
                current = current.get_next()
            return list_
        else :
            return linkedlist

    def list_to_string(self, position):
        pos_in_grid = self.get_pos_in_grid(position)
        self.list_to_string_(pos_in_grid)

    def list_to_string_(self, pos_in_grid):
        pos_list = self.grid[pos_in_grid[0],pos_in_grid[1]]
        if(pos_list != None and pos_list.size()>0):
            current = pos_list.get_head()
            
            print("\nList in position {} : ".format(pos_in_grid))
            while current:
                print(current.get_data().to_string())
                current = current.get_next()

    def to_string(self):
        print("\nParticules in the grid : ")
        for i in range(self.res[0]):
            for j in range(self.res[1]):
                pos_in_grid = [i,j]
                self.list_to_string_(pos_in_grid)
    
# creating the linkedlist we'll use to store our particules inside the grid
    # https://www.codefellows.org/blog/implementing-a-singly-linked-list-in-python/
class Node(object):

    def __init__(self, data=None, next_node=None):
        self.data = data
        self.next_node = next_node

    def get_data(self):
        return self.data

    def get_next(self):
        return self.next_node

    def set_next(self, new_next):
        self.next_node = new_next

class LinkedList(object):
    def __init__(self, head=None):
        self.head = head

    def get_head(self):
        return self.head

    def insert(self, data):
        new_node = Node(data)
        new_node.set_next(self.head)
        self.head = new_node
    
    def size(self):
        current = self.head
        count = 0
        while current:
            count += 1
            current = current.get_next()
        return count

    def search(self, data):
        current = self.head
        found = False
        while current and found is False:
            if current.get_data() == data:
                found = True
            else:
                current = current.get_next()
        if current is None:
            raise ValueError("Data not in list")
        return current


    def delete(self, data):
        current = self.head
        previous = None
        found = False
        while current and found is False:
            if current.get_data() == data:
                found = True
            else:
                previous = current
                current = current.get_next()
        if current is None:
            raise ValueError("Data not in list")
        if previous is None:
            self.head = current.get_next()
        else:
            previous.set_next(current.get_next())

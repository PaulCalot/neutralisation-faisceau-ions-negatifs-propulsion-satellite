# import
import warnings
import numpy as np
from math import pi as math_pi
from dolfin import Point
from random import random, randint
from pprint import pprint
import matplotlib.pyplot as plt
# local imports
from .utils.utils import segment
from .data_structures.vector import MyVector
from .data_structures.grid import Grid

from .utils.integration_schemes import scipy_integrate_solve_ivp, rk4, euler_explicit

class DSMC(object):
    debug = False

    def __init__(self, system, dsmc_params, integration_scheme, f, out_wall = None):
        
        self.dsmc_params = dsmc_params
        self.system = system
        self.grid = system.get_grid()
        self.particles = self.grid.get_all_particles()

        self.integration_scheme = integration_scheme

        self.f = f

        # reacting to collision with walls
        self.walls = self.system.get_walls()
        self.walls_vector = []
        self.init_walls()

        # tracking what happens
        self.collisions_count = 0
        self.collisions_matrix = np.zeros(self.grid.get_grid().shape, dtype = int)
        self.acceptance_rate = 0
        self.total_considered_couples = 0
        self.mean_v_relative = 0
        self.mean_proba = 0

        # temporary in the case we include several types of particles.
        vr_max = self.dsmc_params['vr_max']
        cross_section = np.pi * self.dsmc_params['effective_diameter']*self.dsmc_params['effective_diameter'] # pi d^2
        self.M_ =  (cross_section* vr_max*self.dsmc_params['Ne']*self.dsmc_params['mean_particles_number_per_cell'])/(2*self.dsmc_params['cell_volume'])
        self.vr_max = self.dsmc_params['vr_max']

        # out wall for flux
        self.out_wall = out_wall

    # --------- adding part to the list ----------- #
    def get_number_of_particles(self):
        return len(self.particles)

    def add(self, particle) -> None:
        self.particles.append(particle)

    def get_list_particles(self):
        return self.particles

    def get_scheme(self):
        return self.integration_scheme     
        
    # ------------------------------- Step Function -----------------------------------#

    def step(self, dt, t, f_args):
        debug=False
        if(debug):
            print("Number of parts : {}".format(len(self.particles)))
            plt.figure(figsize=(30,10))
            self.grid.plot()
        list_previous_positions = []
        for i, part in enumerate(self.particles):
            list_previous_positions.append(part.get_2D_pos())
        
        list_count = []
        # First phase DMSC
        list_part_to_delete = []
        for i, part in enumerate(self.particles):
            self._one_particle_step(dt, t, part, f_args)
            #print(part.get_2D_pos(), list_previous_positions[i])
            count = 0
            while(not self.grid.update(part, list_previous_positions[i])):
                count += 1
                # then the particle got out of the box.
                if not self._collision_with_wall(part):
                    # in that case the particle exited the system by the wall
                    # if a print message appears, then this was not suppoed to be the case
                    # in any case : we delete it.
                    list_part_to_delete.append([i,part])
                    break

            if(debug and count == 10):
                print('\n/!\ particle too much reflection /!\ \n')
                print(part.to_string()+" - count : {}".format(count))
                part.plot()
            list_count.append(count)

        if(debug):
            plt.show()
            plt.hist(list_count, bins = 'auto')
            if(len(list_count)>0):
                plt.title("Mean : {} ; std  : {} ; max : {}".format(np.mean(np.array(list_count)),np.std(np.array(list_count)), np.max(np.array(list_count))))
            plt.show()
        # Second phase DMSC
        self.dsmc_collisions(dt)

        # delete all involved particles
        if(debug): print("Number of particles that got out : {}.".format(len(list_part_to_delete)))
        
        for k,particle in enumerate(list_part_to_delete):
            i, part = particle
            # this is a bit messy as the particles list gets smaller and smaller as we delete particles
            # we have to account for that
            self.particles.pop(i-k)
            try:
                self.grid.remove(part)
            except TypeError: 
                self.grid.remove_with_old_pos(part,list_previous_positions[i])
            # this should be enough.

    def _one_particle_step(self, dt, t, part, f_args):
        """ Update particules speed and position when there is no collision. Used in the step method.

        Args:
            dt ([type]): [description]
            t ([type]): [description]
            part_indx ([type]): [description]
            f_args ([type]): [description]
        """

        # particule params
        m=part.get_mass()
        q=part.get_charge()
        pos = part.get_pos() # MyVector object
        speed = part.get_speed()

        f = self.f
        
        args = [m,q]
        f_args = args + f_args

        # new speed and new position
        
        Y=np.array([pos.x, pos.y, pos.z, speed.x, speed.y, speed.z])
        Y_pot = self.integration_scheme(Y, t, dt, f, f_args)
        
        # updating speed and position
        part.set_pos(MyVector(Y_pot[0],Y_pot[1],Y_pot[2]))
        part.set_speed(MyVector(Y_pot[3],Y_pot[4],Y_pot[5]))

    # ----------------------------- Wall collision -------------------------- #
    def _collision_with_wall(self, part):
        debug=False

        pos = part.get_2D_pos()
        speed = part.get_2D_speed()
        radius = part.get_radius()

        min_time, idx_min_time, min_pos_intersect = 1e15 ,-1, None
        wall = None
        liste = []
        for i in range(len(self.walls)):
            t_coll, pos_intersect = self.handler_wall_collision(pos, -1.0*speed, radius, i)
            liste.append((t_coll, pos_intersect, self.walls[i]))

            if(t_coll<min_time):
                idx_min_time = i
                min_time = t_coll
                min_pos_intersect = pos_intersect
                wall = self.walls[i]
        
        # in case the particle went out by the "open wall"
        # another way is just to see if the particle is below / top / to the left / to the right 
        # of the no wall part of the system (we don't include the given wall in that case)
        # and thus we can simply delete the particle if it gets this way
        # but I am not sure it would work better
        if(self.out_wall != None and wall != None and all([l1==l2 for l1, l2 in zip(self.out_wall, wall)])):
            return False

        if(debug): 
            print('Particule : {}'.format(part.to_string()))
            print('Collision time : {}, position : {}'.format(min_time, min_pos_intersect))
            print('with wall {}'.format(wall))
            print('Other walls : ')
            pprint(liste)
        # reflect position / speed of part with respect to wall idx_min_time
        if(min_pos_intersect != None):
            self.reflect_particle(part, min_time, idx_min_time, min_pos_intersect)
            return True
        else:
            #pprint(liste)
            #print("No wall nearby.")
            return False

    def handler_wall_collision(self, position, speed, radius, wall_indx):
        """ Determine if there is a collision between the particule which position, speed and radius 
        are given in parameters and the wall of index wall_indx.
        If there is, it compute the time to collision and update the events table. 
        
        We suppose the particule is caracterized by its position (x,y), its speed (vx, vy) and its radius r.
        The wall is caracterized by its two extremities : (x1,y1), (x2,y2) where x1 <= x2, if x1=x2, then y1 < y2.
        This allows to compute the directing vector of the wall (x2-x1, y2-y1) which has been normalized and 
        stored in *self.walls_vector[wall_indx]*. We note the normalized vector *a*.

        The formula which is used is to compute the possible collision time is : 
            t_coll_1/2 = (-A sgn(B) +/- r)/|B| = (-A +/- r)/B
        where :
            * A = -x sin(theta) + y cos(theta)
            * B = -vx sin(theta) + vy cos(theta)
            * theta = sgn(ay) arccos(ax)
            
        Note that theses times give the moment the disk crosses the infinite line formed from the wall,
        not strictly the wall...

        If B = 0 : we consider that there is not collision and return t_coll = np.nan
        
        If B != O : then a necessary condition is to have both t_coll_1 > 0 and t_coll_2 > 0.
        Indeed : * The first time the particule collides, is when its closest point to the wall collides with it. 
                 * The second time is for when the furthest point to the wall collided with it.
        In such a case, we have to verify that the disk crossing the line occurs on the portion of the line 
        which is the wall. To do that, we compute the position of the particule at the time of collision and verify that it 
        is on the "wall" segment. If it is we return t_coll = min(t_coll_1, t_coll_2). Else, np.nan.

        Args:
            part_indx (int): index of the particule in self.particules
            position (MyVector): position of the particule
            speed (MyVector): speed of the particule
            radius (float): radius of the particule
            wall_indx (int): index of the wall in self.walls

        Returns:
            int, MyVector: the time before the wall and particule collides. Return np.nan is no collision is possible. 
        """
        debug = False
        # p index of the part
        wall_directing_vector = self.walls_vector[wall_indx]
        x1, y1, x2, y2 = self.walls[wall_indx] # x1<=x2 etc.

        # angle
        theta = np.sign(wall_directing_vector.y)*np.arccos(wall_directing_vector.x)

        # A and B
        sTheta, cTheta = np.sin(theta), np.cos(theta)
        B = -speed.x*sTheta+speed.y*cTheta
        if B == 0.0 : 
            if(debug): print('B==0.0')
            return np.nan, None # TODO : should we add a tolerance ? It will never be equals to zero exactly...
        A = -position.x*sTheta+position.y*cTheta
        
        # new position of the wall in the new base
        y1_new_base = -x1*sTheta+y1*cTheta
        if(debug):
            y2_new_base = -x2*sTheta+y2*cTheta
            assert(abs(y1_new_base-y2_new_base)<1e-6)
        A_prime = A-y1_new_base
        # possible collision time :
        t_coll_1 = (-A_prime-radius)/B
        t_coll_2 = (-A_prime+radius)/B
        if(debug): print("Collision time with wall : {} or {}".format(t_coll_1, t_coll_2))
        
        t_intersect = max(t_coll_1, t_coll_2)
        
        if(t_intersect > 0):
            # t_intersect = max(t_coll_1, t_coll_2) # because we are not anticipating them anymore
            # t_intersect = -A_prime/B # the time at which the disk crosses the line if its radius were radius=0.
            pos_intersect = position + t_intersect * speed

            wall_extrimity1_coordinates = MyVector(x1,y1)
            wall_extrimity2_coordinates = MyVector(x2,y2)
            # the reason why were are not using pos_intersect.norm is that it should be a 3D vector.
            dP1, dP2, dP3 = wall_extrimity2_coordinates-wall_extrimity1_coordinates, \
                pos_intersect-wall_extrimity1_coordinates, wall_extrimity2_coordinates-pos_intersect 
            norm_1 = dP1.norm()
            norm_2 = dP2.norm()
            # norm_3 = dP3.norm()
            qty=dP1.inner(dP2)/(norm_1*norm_1) # norm_1 cant be 0 because wall segments are not on same points.
            if(qty < 1 and qty > 0):
                return t_intersect, pos_intersect
            # else:
            #     print("\nQty : {} ".format(qty))
            #     print(self.walls[wall_indx])
            #     print(wall_extrimity1_coordinates)
            #     print(wall_extrimity2_coordinates)
            #     print(pos_intersect)
            #     print('\n')
        else:
            if(debug):
                print(t_intersect)
        # if(t_coll_1 > 0 and t_coll_2 > 0):

        #     t_intersect = -A_prime/B # the time at which the disk crosses the line if its radius were radius=0.
        #     pos_intersect = position + t_intersect * speed 
            
        #     if(self.debug): print("Intersection position with wall : {}, {}".format(pos_intersect.x, pos_intersect.y))

        #     wall_extrimity1_coordinates = MyVector(x1,y1)
        #     wall_extrimity2_coordinates = MyVector(x2,y2)
        #     # the reason why were are not using pos_intersect.norm is that it should be a 3D vector.
        #     dP1, dP2, dP3 = wall_extrimity2_coordinates-wall_extrimity1_coordinates, \
        #         pos_intersect-wall_extrimity1_coordinates, wall_extrimity2_coordinates-pos_intersect 
        #     norm_1 = dP1.norm()
        #     norm_2 = dP2.norm()
        #     # norm_3 = dP3.norm()
        #     if(norm_2 != 0.0):
        #         cosAngle = dP1.inner(dP2)/(norm_1*norm_2)
        #         distance_pos_intersect_to_wall = norm_2*np.sqrt(1-cosAngle*cosAngle)
        #         tol = 1.0 # for now
        #         if(self.debug): print('Distance intersection to wall vs radius : {} vs {}'.format(distance_pos_intersect_to_wall, radius))
        #         if(distance_pos_intersect_to_wall<radius*tol):
        #             if(self.debug):
        #                 print('Next collision in {} s at position {}'.format(min(t_coll_1, t_coll_2),pos_intersect))
        #             return min(t_coll_1, t_coll_2), pos_intersect
        #     else:
        #         # means the two are exactly the same which should virtually never happen
        #         if(self.debug): 
        #             print('Wall extremity one and intersection are the same.')
        #             print('Next collision in {} s at position {}'.format(min(t_coll_1, t_coll_2),pos_intersect))
        #         #return min(t_coll_1, t_coll_2), pos_intersect

            # if(self.debug): print("Wall position : ({}, {}) - ({}, {})".format(x1,y1,x2,y2))
            # if(self.debug): print("Norm : 1, 2, 3: {}, {}, {} \t - \t diff vs radius : {} vs {}.".format(norm_1,norm_2,norm_3,abs(norm_1-norm_2-norm_3), radius))
            # # TODO : should we use a tolerance here too ?
            # if(abs(norm_1-norm_2-norm_3)<radius): # TODO : depending on the computation time - it may not "see" the wall
            #     # because we are too far from the theoretical position ?
            #     # If we want something very stable, we should make the tolerance bigger than radius
            #     if(self.debug):
            #         print('Next collision in {} s at position {}'.format(min(t_coll_1, t_coll_2),pos_intersect))
            #     return min(t_coll_1, t_coll_2), pos_intersect
        if(debug): print('Default out.')
        return np.nan, None

    def reflect_particle(self, part, time, idx, pos_intersect):
        wall_directing_vector = self.walls_vector[idx]
        if(self.debug) : print('Reflected particles')
        # SPEED reflection
        # angle
        theta = float(np.sign(wall_directing_vector.y)*np.arccos(wall_directing_vector.x)) # angle between the directing vector and (0,1)
        # theta must be in degree .... 
        theta = theta*180/np.pi
        # old speed
        old_speed = part.get_2D_speed()
        intermediary_speed = old_speed.rotate(-theta)
        intermediary_speed_2 = MyVector(intermediary_speed.x, -intermediary_speed.y)
        new_speed_2d = intermediary_speed_2.rotate(theta)
        
        # setting new speed 
        new_speed_3d = MyVector(new_speed_2d.x, new_speed_2d.y, part.get_speed().z)
        part.set_speed(new_speed_3d)

        # POSITION reflection
        old_pos = part.get_pos()
        new_pos = MyVector(pos_intersect.x, pos_intersect.y, part.get_pos().z) + (0.1*(0.5-random())+1)*time*new_speed_3d
        part.set_pos(new_pos)

        if(self.debug):
            print("Old position : {}".format(old_pos))
            print("New position : {}".format(new_pos))

    # ---------------------------- Collision phase -------------------------- #
    def dsmc_collisions(self, dt) -> None:
        grid = self.grid.get_grid()
        M = self.M_*dt
        for idx_line, line in enumerate(grid):
            for idx_cell, cell in enumerate(line):
                if(cell == None):
                    continue
                
                # may i should do differently # and (self.sparsed_space[i][j] == 1 if self.use_sparsed_space else True)
                nc = cell.get_size()
                mcand = int(nc*M) # number of pairs to select to compute
                self.total_considered_couples += mcand
                
                self.handle_one_cell_collisions(mcand, nc, cell, idx_line, idx_cell)

    def handle_one_cell_collisions(self, Mcand, Nc, cell, idx_line, idx_cell):
        for k in range(Mcand):
            i,j = randint(0,Nc-1), randint(0,Nc-1)
            while(i==j):
                j = randint(0,Nc-1)

            # at this point we don't assert that we are not computing several time the same pairs
            # but since Mcand << Nc in theory, we should be fine most of the time

            part1, part2 = cell.get_pair(i, j) # this is O(N) which is pretty long
            v_r_norm = (part2.get_speed()-part1.get_speed()).norm()

            #If this relative speed norm is higher than the one we set so far, we update.
            if(v_r_norm>self.vr_max): # remember for next time
                self.vr_max = v_r_norm

            r = random() # uniform deviate in (0,1)
            acceptance_proba = v_r_norm/self.vr_max

            if(r<acceptance_proba):
                self.update_speed_dsmc(part1, part2, v_r_norm)

                # tracking
                self.collisions_matrix[idx_line,idx_cell] += 1
                self.acceptance_rate += 1
                self.collisions_count+=1
            # tracking
            self.mean_proba += acceptance_proba
            self.mean_v_relative += v_r_norm


    def update_speed_dsmc(self, part1, part2, v_r_norm):
        r = random() # uniform deviate in (0,1)
        # for now, no loss of energy
        # page 5/7 of the Direct Simulation MC method paper
        q = 2*r-1
        cTheta = q
        sTheta = np.sqrt(1-q*q)
        phi = 2*math_pi*r

        # computations of the speed
        v_cm = 0.5*(part1.get_speed()+part2.get_speed()) # this quantity is conserved
        v_r_ = v_r_norm*(MyVector(sTheta*np.cos(phi), sTheta*np.sin(phi), cTheta))
        
        # setting the new speeds
        part1.set_speed(v_cm+0.5*v_r_) 
        part2.set_speed(v_cm-0.5*v_r_)

    # ------------------------------- Init ---------------------------- #

    def add_particles(self, particles_list):
        tmp = 0
        for part in particles_list:
            part.set_id(tmp)
            tmp+=1
            init_position = part.get_2D_pos()
            init_speed = part.get_2D_speed()
            if(self.debug): print("\n\n"+part.to_string())
            B = True
            while(not self.grid.add_and_verify(part)): # adding it to the grid
                # then the particle has been initialized out of the box.
                if not self._collision_with_wall(part):
                    # in that case the particle exited the system by the wall
                    # if a print message appears, then this was not suppoed to be the case
                    # in any case : we delete it.
                    print('ERROR - newly initialized particle can not be linked to any wall it could have exited')
                    print(part.to_string())
                    print('Initial position : {}'.format(init_position))
                    print('Init speed : {}'.format(init_speed))
                    part.plot()
                    B = False
                    break 
            if(B):             
                self.add(part)
                
    def init_walls(self):
        for i, wall in enumerate(self.walls):
            x1,y1,x2,y2 = wall
            min_x, max_x, min_y, max_y = x1, x2, y1, y2
        
            if(x1 > x2):
                min_x, max_x = x2, x1
                min_y, max_y = y2, y1
            elif(x1 < x2):
                min_x, max_x = x1, x2
                min_y, max_y = y1, y2
            else :
                if(y1 > y2):
                    min_x, max_x = x2, x1
                    min_y, max_y = y2, y1
                elif(y2 > y1):
                    min_x, max_x = x1, x2
                    min_y, max_y = y1, y2
                else :
                    warnings.warn("The {}th segment-wall has the same two extremities : (x1,y1) = ({},{}) and (x2,y2) = ({},{}).".format(i,x1,y1,x2,y2))

            #self.walls[i] = segment(Point(min_x, min_y),Point(max_x, max_y)) # not updating walls anymore as it was creating issues...
            self.walls_vector.append(MyVector(max_x-min_x, max_y-min_y).normalize())# the directing vectors

    # ------------------------------- Getter and setter -------------------------- #
   
    def get_collisions_count(self):
        return self.collisions_count
    
    def get_acceptance_rate(self):
        return self.acceptance_rate/self.total_considered_couples

    def get_mean_vr_norm(self):
        return self.mean_v_relative/self.total_considered_couples

    def get_vr_norm(self):
        # return self.dsmc_params['vr_max']
        return self.vr_max
    
    def get_mean_proba(self):
        return self.mean_proba/self.total_considered_couples
    
    # ------------------ Saving function ------------------ # 
    # Not sure I should put it there ...

    def save_collisions_matrix(self, name, iteration = 0):
        from os.path import isfile
        if(isfile(name)):
            mode = 'a'
        else:
            mode = 'w'
        with open(name, mode = mode) as txt_file:
            head = "Iteration {} :\n".format(str(iteration))
            txt_file.write(head)
            L=[]
            #print(self.collisions_matrix)
            for k in range(self.collisions_matrix.shape[0]):
                line = self. collisions_matrix[k]
                list_ = [str(line[nb]) for nb in range(len(line))] + ['\n']
                L.insert(0, '  '.join(list_))
            for string in L:
                txt_file.write(string)
            txt_file.write('\n')


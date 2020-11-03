# import
import warnings
import numpy as np

from dolfin import Point

# local imports
from .utils import segment
from .vector import MyVector

class CollisionHandler(object):
    debug = False
    def __init__(self, particules, walls, f, eta = 0, p = 0):
        # TODO : eta and p are for lose of speed or lose of charge when contact with wall or particules
        self.particules = particules
        self.walls = walls
        self.nb_parts = len(particules)
        self.nb_walls = len(walls)
        self.walls_vector = []
        # update speed function
        self.f = f 
        # eta : coefficient that caracterize the energy loss - eta = 0 : no loss
        self.eta = eta
        # probability of losing a charge for the particule when colliding
        self.lose_charge_proba = p 

        # and sorting segments point by incrementing x
        # /!\ useful ? Pretty sure it's not.
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
            self.walls[i] = segment(Point(min_x, min_y),Point(max_x, max_y))

            self.walls_vector.append(MyVector(max_x-min_x, max_y-min_y).normalize())# the directing vectors
        
        #self.walls_vector = np.array(self.walls_vector, dtype = MyVector)

        # events table
        self.events = np.zeros((self.nb_parts, self.nb_parts + self.nb_walls), dtype = float) # in the columns, walls come first
        for i in range(self.nb_parts): self.update_events(i)

    # ------------------- Updating events functions ------------------------- #

    def update_all_events(self):
        for i in range(self.nb_parts): self.update_events(i)

    def update_events(self, p):
        """ Update the events of the p-th disk.
        Note that *self.events* array is a 2D-matrix that has *nb_walls+nb_parts* columns (where the first nb_walls ones are for walls) and
        *nb_parts* lines.

        Args:
            p (int): index of the particule in self.particules.
        """
        part = self.particules[p]

        # current particule pos and speed :
        pos = part.get_pos()
        speed = part.get_speed()

        # radius
        r = part.get_radius()

        # no self-collision !
        self.events[p][p+self.nb_walls] = np.nan
        
        # collisions with a wall
        for i in range(self.nb_walls):
            self.events[p,i] = self.handler_wall_collision(pos, speed, r, i)

        # collisions with other disks
        for i in range(0, self.nb_parts):
            # updating only if particule i is not p.
            if i != p:
                t = self.handler_particules_collision(p, i)
                self.events[p,i+self.nb_walls] = t
                self.events[i,p+self.nb_walls] = t
    
        if(self.debug): print("\nUpdated events of particule {} : {}\n".format(p,self.events))


    
    # ------------------------------- Step Function -----------------------------------#
    def step(self, dt, t, f_args):
        #if(self.debug): print(self.events)
        try :
            ind = np.unravel_index(np.nanargmin(self.events), self.events.shape)
            t_coll = self.events[ind]
        except :
            ind = 0
            t_coll = np.nan
        if(self.debug): print("Current time : {} sec. Time before next collision : {} sec.".format(t,t_coll))
        if dt < t_coll or np.isnan(t_coll):
            for part_indx in range(len(self.particules)):
                self._one_particule_step(dt, t, part_indx, f_args)
            self.events -= dt
        else:
            # TODO : make sure there is no problem here with updating the position and speed of the particule WHEREAS there will be a collision.
            # i.e. We certainly don't want to have the particule going to far in the wall as it will "keep bouncing" in such a case.
            for part_indx in range(len(self.particules)):
                self._one_particule_step(t_coll, t, part_indx, f_args)
            self.events -= t_coll
            self._new_velocities_collision(ind) # update velocity of part of index ind
            self.step(dt - t_coll, t+t_coll, f_args)
    
    # ----------------------------- Collision handlers ------------------------------- #

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
        # p index of the part
        wall_directing_vector = self.walls_vector[wall_indx]
        x1, y1, x2, y2 = self.walls[wall_indx] # x1<=x2 etc.

        # angle
        theta = np.sign(wall_directing_vector.y)*np.arccos(wall_directing_vector.x)

        # A and B
        sTheta, cTheta = np.sin(theta), np.cos(theta)
        B = -speed.x*sTheta+speed.y*cTheta
        if B == 0.0 : return np.nan# TODO : should we add a tolerance ? It will never be equals to zero exactly...
        A = -position.x*sTheta+position.y*cTheta
        
        # new position of the wall in the new base
        y1_new_base = -x1*sTheta+y1*cTheta
        if(self.debug):
            y2_new_base = -x2*sTheta+y2*cTheta
            assert(abs(y1_new_base-y2_new_base)<1e-6)
        A_prime = A-y1_new_base
        # possible collision time :
        t_coll_1 = (-A_prime-radius)/B
        t_coll_2 = (-A_prime+radius)/B
        if(self.debug): print("\nCollision time with wall : {} or {}".format(t_coll_1, t_coll_2))

        if(t_coll_1 > 0 and t_coll_2 > 0):
            t_intersect = -A_prime/B # the time at which the disk crosses the line if its radius were radius=0.
            pos_intersect = position + t_intersect * speed 
            
            if(self.debug): print("Intersection position with wall : {}, {}\n".format(pos_intersect.x, pos_intersect.y))

            wall_extrimity1_coordinates = MyVector(x1,y1)
            wall_extrimity2_coordinates = MyVector(x2,y2)
            # the reason why were are not using pos_intersect.norm is that it should be a 3D vector.
            dP1, dP2, dP3 = wall_extrimity2_coordinates-wall_extrimity1_coordinates, \
                pos_intersect-wall_extrimity1_coordinates, wall_extrimity2_coordinates-pos_intersect
            norm_1 = dP1.norm()
            norm_2 = dP2.norm()
            norm_3 = dP3.norm()
            if(self.debug): print("Wall position : ({}, {}) - ({}, {})".format(x1,y1,x2,y2))
            if(self.debug): print("Norm : 1, 2, 3 : {}, {}, {}\n".format(norm_1,norm_2,norm_3))
            # TODO : should we use a tolerance here too ?
            if(abs(norm_1-norm_2-norm_3)<radius): # TODO : depending on the computation time - it may not "see" the wall
                # because we are too far from the theoretical position ?
                # If we want something very stable, we should make the tolerance bigger than radius
                return min(t_coll_1, t_coll_2)
        return np.nan
    
    # TODO : it could be best to add these functions to an external class so we can choose how we handle collision
    # independently from the algorithm that handles them (with the events table)
    
    def handler_particules_collision(self, p, i):
        """Helper method to get the collision time between two disks of radius r"""

        # TODO : make it work with disks of different radius ? I think that we can approximate I, I-, I+  to the same radius.
        
        # getting particules positions and speed
        part1 = self.particules[p]
        part2 = self.particules[i]

        dr = part2.get_2D_pos() - part1.get_2D_pos()
        dv = part2.get_2D_speed() - part1.get_2D_speed()
        if(dv.x == 0 and dv.y == 0.0):
            return np.nan
        dvdr = np.dot(dv,dr)
        
        radius = part1.get_radius()

        # the disks need to approach each other to collide
        if dvdr > 0:
            return np.nan
        
        # check if they will collide
        sigma = 4*radius**2
        dvdv = np.dot(dv,dv)
        drdr = np.dot(dr,dr)
        beta = 1-dvdv/(dvdr*dvdr)*(drdr-sigma)
        if beta<0:
            return np.nan
        
        # return time to collision
        alpha = dvdr/dvdv
        return -alpha*(1-np.sqrt(beta))
    
    # ------------------------------------ Velocity and/or position updates --------------------------------#
    def _one_particule_step(self, dt, t, part_indx, f_args):
        """ Update particules speed and position when there is no collision. Used in the step method.

        Args:
            dt ([type]): [description]
            t ([type]): [description]
            part_indx ([type]): [description]
            f_args ([type]): [description]
        """
        particule=self.particules[part_indx]

        # particule params
        m=particule.get_mass()
        q=particule.get_charge()
        espece=particule.get_part_type()
        pos = particule.get_pos() # MyVector object
        speed = particule.get_speed()

        f = self.f
        # new speed and new position
        Y=np.array([pos.x, pos.y, pos.z, speed.x, speed.y, speed.z])
        k1=np.array(f(Y,t,m,q,*f_args))
        k2=np.array(f(Y+.5*dt*k1, t+.5*dt,m,q,*f_args))
        k3=np.array(f(Y+.5*dt*k2, t+.5*dt,m,q,*f_args))
        k4=np.array(f(Y+dt*k3, t+dt,m,q,*f_args))
        Y_pot=Y+dt*(1/6*k1+1/3*k2+1/3*k3+1/6*k4)

        # updating speed and position
        particule.set_pos(MyVector(Y_pot[0],Y_pot[1],Y_pot[2]))
        particule.set_speed(MyVector(Y_pot[3],Y_pot[4],Y_pot[5]))


    def _new_velocities_collision(self, ind):
        """ Update velocity in case of a collision. Does not update speed however.

        Args:
            ind ([type]): [description]
        """
        if(self.debug): print("Particule {} is colliding with wall {}".format(ind[0],ind[1]) if ind[1]<self.nb_walls else \
            "Particule {} is colliding with particule {}".format(ind[0],ind[1]-self.nb_walls))

        # TODO : make the collision more "realistic" : 3D collision, taking into account the angle for inelastic collision + best way to do it?
        # ind is the position in event matrix where the next collision occures
        p1 = ind[0] # the particule indx in the matrix
        part1 = self.particules[p1]
        
       
        # Event was a collision with wall
        if ind[1] < self.nb_walls:
            wall_directing_vector = self.walls_vector[ind[1]]

            # angle
            theta = float(np.sign(wall_directing_vector.y)*np.arccos(wall_directing_vector.x)) # angle between the directing vector and (0,1)
            # theta must be in degree .... 
            theta = theta*180/np.pi
            # old speed 
            old_speed = part1.get_2D_speed()
            intermediary_speed = old_speed.rotate(-theta)
            intermediary_speed_2 = MyVector(intermediary_speed.x, -intermediary_speed.y)
            new_speed = intermediary_speed_2.rotate(theta)
            
            # setting new speed
            part1.set_2D_speed(new_speed)

            part1.set_2D_speed(part1.get_2D_speed()*(1-self.eta))

            if(np.random.random_sample()<self.lose_charge_proba):
                part1.lose_charge()
            
            if(self.debug): print("New SPEED : [{}] - pos : {} - speed : {}".format(p1,(round(part1.get_2D_pos().x,3),round(part1.get_2D_pos().y,3)),\
                (round(part1.get_2D_speed().x,3),round(part1.get_2D_speed().y,3))))

            self.update_events(p1) # same here as for the wall...
            
        # Event was collision with another disk
        else:
            
            p2 = ind[1]-self.nb_walls
            part2 = self.particules[p2]
            sigma = 2 * part1.get_radius()

            # old speeds
            old_speed1 = part1.get_2D_speed()
            old_speed2 = part2.get_2D_speed()

            delta_v = part2.get_2D_speed() - part1.get_2D_speed()
            delta_r = part2.get_2D_pos() - part1.get_2D_pos()

            dvdr = np.dot(delta_v, delta_r)

            imp = 2*1*1*dvdr/(sigma*(1+1))
            imp_v = imp*delta_r/sigma
        
            part1.set_2D_speed(part1.get_2D_speed() + imp_v)
            part2.set_2D_speed(part2.get_2D_speed() - imp_v)

            # in the case of a two-disk collisions, we should have something different
            # as for the charge loss etc.
            part1.set_2D_speed(part1.get_2D_speed()*(1-self.eta))

            if(np.random.random_sample()<self.lose_charge_proba):
                part1.lose_charge()

            part2.set_2D_speed(part2.get_2D_speed()*(1-self.eta))

            if(np.random.random_sample()<self.lose_charge_proba):
                part2.lose_charge()

            self.update_events(p1) # same here as for the wall...
            self.update_events(p2)
            if(self.debug): 
                print("New SPEED : [{}] - pos : {} - speed : {}".format(p1,(round(part1.get_2D_pos().x,3),round(part1.get_2D_pos().y,3)),\
                (round(part1.get_2D_speed().x,3),round(part1.get_2D_speed().y,3))))
                print("New SPEED : [{}] - pos : {} - speed : {}".format(p2,(round(part2.get_2D_pos().x,3),round(part2.get_2D_pos().y,3)),\
                (round(part2.get_2D_speed().x,3),round(part2.get_2D_speed().y,3))))
    
    # ------------------------------- Getter and setter -------------------------- #

    def get_next_collision(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        try : 
            ind = np.unravel_index(np.nanargmin(self.events), self.events.shape)
            t_coll = self.events[ind]
        except :
            t_coll = np.nan 
        return t_coll
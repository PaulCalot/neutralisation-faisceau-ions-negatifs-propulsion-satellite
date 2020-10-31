# import
import warnings
import numpy as np

from dolfin import Point
from .utils import segment
# local imports
from .vector import MyVector


class CollisionHandler(object):
    
    tolerance = 1e-12 # depends on the size of the atom

    def __init__(self, particules, walls, f, eta, p):
        # TODO : eta and p are for lose of speed or lose of charge when contact with wall or particules
        self.particules = particules
        self.walls = walls
        self.nb_parts = len(particules)
        self.nb_walls = len(walls)
        self.walls_coeff = np.zeros((self.nb_walls,2))
        self.f = f 
        self.eta = eta
        self.lose_charge_proba = p
        # precomputing coefficients for the walls
        # and sorting segments point by incrementing x
        for i, wall in enumerate(walls):
            x1,y1,x2,y2 = wall
            if(x1 == x2):
                a = 0
                b = x1
            else :
                a = (y2-y1)/(x2-x1) if x2!=x1 else 0
                b = y1 - a * x1
            
            self.walls_coeff[i][0] = a
            self.walls_coeff[i][1] = b

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

        # events table
        self.events = np.array((self.nb_parts, self.nb_parts + self.nb_walls), dtype = float) # in the columns, walls come first

    def update_events(self, p):
        """Update the events for the p^th disk"""
    
        # getting particules positions and speed
        pos = [part.get_2D_pos() for part in self.particules]
        speed = [part.get_2D_speed() for part in self.particules]

        # current particule pos and speed :
        x, y, vx, vy = pos[p].x, pos[p].y, speed[p].x, speed[p].y
        r = self.particules[p].get_radius()

        # no self-collision!
        self.events[p,p+self.nb_walls] = np.nan
        
        # collisions with a wall
        for i in range(self.nb_walls):

            # make a function that handles the collision with a wall
            a, b = self.walls_coeff[i]
            A = vy + a * vx
            B = y - a * x - b
            t_coll1 = (B + r)/A
            t_coll2 = (B - r)/A
            if (t_coll1>0 and t_coll2>0):
                t_coll = min(t_coll1, t_coll2)
                t_intersect = B/A
                # we still have to verify that we indeed have a solution
                # (weather or not we are on the segment indeed)
                x_sol, y_sol = x+t_intersect*vx, y+t_intersect*vy
                x1,y1,x2,y2 = self.walls[i]

                if((abs(x1)-self.tolerance < x_sol and abs(x2)+self.tolerance > x_sol) or \
                    (abs(y1)-self.tolerance < y_sol and abs(y2)+self.tolerance > y_sol)):
                    self.event[p,i] = t_coll
                else :
                    self.event[p,i] = np.nan
            
            else : 
                t_coll = max(t_coll1, t_coll2) 
                # this is tricky : basically, if we find ourselves in such a case, 
                # it's most likely because we missed the collision (which is -min(t_coll1, t_coll2) "behind" us)
                # in this case, the right thing is most likely to compute the collision right away ? 
                # I think it just should not append however unless a particule starts too close from a wall.
                # that's why I set "0" in this case.
                if(t_coll >= 0):
                    self.event[p,i] = 0
                else : 
                    self.event[p,i] = np.nan
    
        # collisions with other disks
        for i in range(self.nb_walls, self.n_disks + self.nb_walls):
            if i != p:
                t = self._time_to_collision(p, i)
                self.events[p,i+self.nb_walls] = t
                self.events[i,p+self.nb_walls] = t
    


    def update_particule(self, dt, t, part_indx, E, zone):
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
        k1=np.array(f(Y,t,m,q,zone,E))
        k2=np.array(f(Y+.5*dt*k1, t+.5*dt,m,q,zone,E))
        k3=np.array(f(Y+.5*dt*k2, t+.5*dt,m,q,zone,E))
        k4=np.array(f(Y+dt*k3, t+dt,m,q,zone,E))
        Y_pot=Y+dt*(1/6*k1+1/3*k2+1/3*k3+1/6*k4)

        # updating speed and position
        particule.set_pos(MyVector(Y_pot[0],Y_pot[1],Y_pot[2]))
        particule.set_speed(MyVector(Y_pot[3],Y_pot[4],Y_pot[5]))

    def step(self, dt, t, E, zone):
        ind = np.unravel_index(np.nanargmin(self.events), self.events.shape)
        t_coll = self.events[ind]
        
        if dt < t_coll:
            for part_indx in range(len(self.particules)):
                self.update_particule(dt, t, part_indx, E, zone)
                self.events -= dt
        else:
            for part_indx in range(len(self.particules)):
                self.update_particule(t_coll, t, part_indx, E, zone)
                self.events -= t_coll
            self._new_velocities(ind) # update velocity of part of index ind
            self.step(dt - t_coll)
    
    # TODO : it could be best to add these functions to an external class so we can choose how we handle collision
    # independently from the algorithm that handles them (with the events table)
    
    def _time_to_collision(self, p, i):
        """Helper method to get the collision time between two disks of radius r"""
        # TODO : make it work with disks of different radius ? I think that we can approximate I, I-, I+  to the same radius.
        # TODO : get_pos -> get_pos_2D
        # getting particules positions and speed
        part1 = self.particules[p]
        part2 = self.particules[i]

        dr = part2.get_2D_pos() - part1.get_2D_pos()
        dv = part2.get_2D_speed() - part1.get_2D_speed()
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
    
    def _new_velocities(self, ind):

        # TODO : make the collision more "realistic" : 3D collision, taking into account the angle for inelastic collision + best way to do it?
        # ind is the position in event matrix where the next collision occures
        p1 = ind[0] # the particule indx in the matrix
        part1 = self.particules[p1]
        
        # Event was a collision with wall
        if ind[1] < self.nb_walls:
            a, b = self.walls_coeff[ind[1]]
            direction_vector = MyVector(a,1).normalize()
            normal = MyVector(-a,1).normalize()
            old_speed = part1.get_2D_speed()
            theta = 2*np.arcos(direction_vector.inner(old_speed.normalize()))
            dot_product =  normal.inner(old_speed.normalize())
            if(dot_product>0):
                theta = - theta
            part1.rotate_speed_2D(theta) # rotate speed along +z axis. 
            
            # at least we have theta so this should be ok if we want to take it into acount
            # not sure about this one (we lose speed weird ...)
            part1.set_2D_speed(part1.get_2D_speed()*(1-self.eta))

            if(np.random.random_sample()<self.lose_charge_proba):
                part1.lose_charge()

            self.update_events(p1)

        # Event was collision with another disk
        else:
            
            p2 = ind[1]- self.nb_walls
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

            self.update_events(p1)
            self.update_events(p2)
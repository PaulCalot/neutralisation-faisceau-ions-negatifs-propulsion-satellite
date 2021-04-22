import numpy as np
cimport numpy as np

DTYPE = np.short

cimport cython
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def handler_wall_collision(double px, double py, double vx, double vy, double radius, double x1, double x2, double y1, double y2, double ax, double ay):
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
        cdef double sTheta, cTheta
        cdef double A_prime, B
        # cdef double y1_new_base
        # cdef double t_coll_1, t_coll_2
        cdef double t_intersect
        cdef double pix, piy
        cdef double dP1x, dP1y, dP2x, dP2y
        cdef double qty
        cdef double t1, t2

        # A and B
        if(x2-x1 > 0):
            sTheta, cTheta = y2-y1, x2-x1
        else:
            sTheta, cTheta = y1-y2, x2-x1

        B = -vx*sTheta+vy*cTheta
        if B == 0.0 : 
            return 1e15, None, None # TODO : should we add a tolerance ? It will never be equals to zero exactly...
        
        # A = -position[0]*sTheta+position[1]*cTheta
        # new position of the wall in the new base
        # y1_new_base = -wall[0]*sTheta+wall[1]*cTheta
        A_prime =  -px*sTheta+py*cTheta+x1*sTheta-y1*cTheta
        # possible collision time :
        t1 = (-A_prime-radius)/B
        t2 = (-A_prime+radius)/B
        if(t1>t2):
            t_intersect = t1
        else:
            t_intersect = t2

        
        if(t_intersect > 0):
            # t_intersect = max(t_coll_1, t_coll_2) # because we are not anticipating them anymore
            # t_intersect = -A_prime/B # the time at which the disk crosses the line if its radius were radius=0.
            pix = px + t_intersect * vx
            piy = py + t_intersect * vy

            # the reason why were are not using pos_intersect.norm is that it should be a 3D vector.
            dP1x = x2-x1
            dP1y = y2-y1
            dP2x = x1-pix
            dP2y = y1-piy

            qty=(dP1x*dP2x+dP1y*dP2y)/(dP1x*dP1x+dP1y*dP1y) # norm_1 cant be 0 because wall segments are not on same points. # np.vdot(dP1,dP2)
            if(qty < 1 and qty > 0):
                return t_intersect, pix, piy
        return 1e15, None, None
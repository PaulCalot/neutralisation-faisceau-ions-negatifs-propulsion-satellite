import numpy as np
cimport numpy as np

DTYPE = np.short

cimport cython
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def handler_wall_collision(double px, double py, double vx, double vy, double radius, double x1, double y1, double x2, double y2, double ax, double ay):

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
            #print("B=0")
            return 1e15, None, None # TODO : should we add a tolerance ? It will never be equals to zero exactly...
        
        # A = -position[0]*sTheta+position[1]*cTheta
        # new position of the wall in the new base
        # y1_new_base = -wall[0]*sTheta+wall[1]*cTheta
        A_prime =  -px*sTheta + py*cTheta + x1*sTheta - y1*cTheta
        #print(f"A_prime = {A_prime}")
        # possible collision time :
        t1 = (-A_prime-2*radius)/B
        t2 = (-A_prime+2*radius)/B

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
            dP2x = pix-x1
            dP2y = piy-y1

            qty=(dP1x*dP2x+dP1y*dP2y)/(dP1x*dP1x+dP1y*dP1y) # norm_1 cant be 0 because wall segments are not on same points. # np.vdot(dP1,dP2)
            if(qty < 1 and qty > 0):
                return t_intersect, pix, piy
        #print("Default")
        return 1e15, None, None
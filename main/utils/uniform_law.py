# -*- coding: utf-8 -*-
"""
Spyder Editor

Link :
    https://en.wikipedia.org/wiki/Spherical_cap#cite_note-handbook-1
    Problem I have : https://math.stackexchange.com/questions/936681/volume-of-a-sphere-corner
    
"""
from math import pi
import numpy as np
import matplotlib.pyplot as plt

def pdf(v,a):
    if(v<0):
        return 0.0
    elif(v < np.sqrt(3.0)*a/2.0):
        sphere_volume = 4/3.0*pi*v*v*v
        cube_volume = a*a*a
        m = a/2.0
        if(v<m):
            return sphere_volume/cube_volume
        else:
            r, h = v, v-m
            spherical_cap_volume = pi*h*h/3.0*(3*r-h)    
            if(v<np.sqrt(2.0)*m):
                return (sphere_volume-6*spherical_cap_volume)/cube_volume
            else:
                # at this point, the spherical caps start to overlap...
                mp = np.sqrt(h*(2*r-h))
                L = mp-m
                v_corr_1 = L*L*0.5*(2*np.sqrt(mp*mp-m*m))  # still not exactly that though 
                # Pb : la partie en plus ne peut se décomposer avec un cylindre ...
                v_corr = v_corr_1
                return (sphere_volume-6*spherical_cap_volume+12*v_corr)/cube_volume
    else :
        return 1
    
a=2.0
V = np.linspace(-0.1,2,10000)
pdf_speed_norm = np.array([pdf(v,a) for v in V])
plt.plot(V,pdf_speed_norm)
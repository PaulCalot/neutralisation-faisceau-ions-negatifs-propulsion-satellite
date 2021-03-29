from __future__ import print_function

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import dolfin
import mshr
import numpy as np
import scipy.integrate as integrate
from fenics import *
# TODO : delete unused imporrt


def segment(Point1,Point2):
    """
    Renvoie une liste des 4 coordonnées des extrémités à partir des deux points donnés
    [x1,y1,x2,y2]
    """
    return [Point1[0],Point1[1],Point2[0],Point2[1]]


def get_mesh(mesh_dict):
    """Returns the mesh and the segment grid.

    Args:
        L_mot (int): 
        l_mot ([type]): [description]
        h_grid ([type]): [description]
        L_vacuum ([type]): [description]
        l_vacuum ([type]): [description]
        L_1 ([type]): [description]
        l_1 ([type]): [description]
        L_2 ([type]): [description]
        l_2 ([type]): [description]
        L_3 ([type]): [description]
        l_3 ([type]): [description]
        delta_vert_12 ([type]): [description]
        delta_vert_23 ([type]): [description]
        mesh_resolution ([type]): [description]
        refine_mesh ([type]): [description]

    Returns:
        [type]: [description]
    """
    L_mot = mesh_dict['L_mot']
    l_mot = mesh_dict['l_mot']
    L_vacuum = mesh_dict['L_vacuum']
    l_vacuum = mesh_dict['l_vacuum']
    L_1 = mesh_dict['L_1']
    l_1 = mesh_dict['l_1']
    L_2 = mesh_dict['L_2']
    l_2 = mesh_dict['l_2']
    L_3 = mesh_dict['L_3']
    l_3 = mesh_dict['l_3']
    delta_vert_12 = mesh_dict['delta_vert_12']
    delta_vert_23 = mesh_dict['delta_vert_23']
    mesh_resolution = mesh_dict['mesh_resolution']
    refine_mesh = mesh_dict['refine_mesh']

    h_grid = l_1 + l_2 + l_3 + delta_vert_12 + delta_vert_23

    P1 = Point(-L_mot/2, -l_mot/2)
    P2 = Point(L_mot/2, l_mot/2)
    P3 = P1 + Point(L_mot, -h_grid)
    P4 = P1 + Point(-(L_vacuum-L_mot)/2, -(h_grid+l_vacuum))
    P5 = P3 + Point((L_vacuum-L_mot)/2, 0)

    moteur = mshr.Rectangle(P1, P2)
    grid = mshr.Rectangle(P1, P3)
    vacuum = mshr.Rectangle(P4, P5)

    zone = moteur + grid + vacuum

    P6 = P1 + Point(L_1/2, -l_1)
    P7 = P6 + Point(-L_1/2 , -delta_vert_12)
    P8 = P7 + Point(L_2/2, -l_2)
    P9 = P8 + Point(-L_2/2 , -delta_vert_23)
    P10 = P9 + Point(L_3/2, -l_3)
    P11 = P3 + Point(-L_3/2, l_3)
    P12 = P11 + Point(L_3/2, delta_vert_23)
    P13 = P12 + Point(-L_2/2, l_2)
    P14 = P13 + Point(L_2/2, delta_vert_12)
    P15 = P14 + Point(-L_1/2, l_1)

    rect_1 = mshr.Rectangle(P1, P6)
    rect_2 = mshr.Rectangle(P7, P8)
    rect_3 = mshr.Rectangle(P9, P10)
    rect_4 = mshr.Rectangle(P3, P11)
    rect_5 = mshr.Rectangle(P12, P13)
    rect_6 = mshr.Rectangle(P14, P15)

    zone -= (rect_1 + rect_2 + rect_3 + rect_4 + rect_5 + rect_6)

    segments_list=[
        segment(P1,P2+Point(-L_mot,0)), #0
        segment(P2+Point(-L_mot,0),P2), 
        segment(P2, P15+Point(L_1/2,0)), #2
        segment(P15+Point(L_1/2,0), P15),
        segment(P15, P15+Point(0,-l_1)),
        segment(P15+Point(0,-l_1), P14),
        segment(P14, P14+Point(0,-delta_vert_12)), #6
        segment(P14+Point(0,-delta_vert_12), P13),
        segment(P13, P13+Point(0, -l_2)),
        segment(P13+Point(0, -l_2), P12),
        segment(P12, P12+Point(0,-delta_vert_23)), #10
        segment(P12+Point(0,-delta_vert_23), P11),
        segment(P11, P11+Point(0, -l_3)),
        segment(P11+Point(0, -l_3), P5),
        segment(P5, P5+Point(0,-l_vacuum)),
        segment(P5+Point(0,-l_vacuum), P4),
        segment(P4, P4+Point(0,l_vacuum)),
        segment(P4+Point(0,l_vacuum), P10),
        segment(P10, P10+Point(0, l_3)),
        segment(P10+Point(0, l_3), P9),
        segment(P9, P9+Point(0,delta_vert_23)),#20
        segment(P9+Point(0,delta_vert_23), P8), 
        segment(P8, P8+Point(0, l_2)),
        segment(P8+Point(0, l_2), P7),
        segment(P7, P7+Point(0,delta_vert_12)),#24
        segment(P7+Point(0,delta_vert_12), P6), 
        segment(P6, P6+Point(0,l_1)),
        segment(P6+Point(0,l_1),P1)
    ]
    
    mesh=mshr.generate_mesh(zone, mesh_resolution)
    

    
    if(refine_mesh):
        d = mesh.topology().dim()
        
        class To_refine(SubDomain):
            def inside(self, x, on_boundary):
                return x[1]<=0 and x[1]>= -l_mot/2-h_grid-l_vacuum/4

        to_refine = To_refine()
        marker = MeshFunction("bool", mesh, d, False)
        to_refine.mark(marker, True)
        mesh = refine(mesh,marker)

    return mesh, segments_list, zone


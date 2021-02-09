from main import get_mass_part, Particule, MyVector
from main import available_particles
from main import get_gaussian_params_maxwellian, available_particles
# imports,
from dolfin import Point
from random import gauss, random
import numpy as np
from scipy.stats import maxwell, norm
from pprint import pprint

def init_particles_in_system(particles_types, particles_densities, particles_mean_number_per_cell, speed_type, speed_param1,\
     speed_param2, space_size, space_resolution, zone = None, offsets = [0,0], verbose = False, debug = False, *args, **kwargs):
    list_particles=[]
    e = 1.6e-19
    mean = 0
    number_ = 1
    for k in range(len(particles_densities)):
        if(debug):
            print(particles_types)
        type_ = particles_types[k]
        number_ = int(particles_mean_number_per_cell[k]*space_resolution[0]*space_resolution[1])
        speed_init_type = speed_type[k] # uniform or maxwellian
        m, M = speed_param1[k], speed_param2[k]
        if(debug): 
            print('Creating {} particles of type {}.'.format(number_, type_))
            print('Speed type init is {} with params {}, {}.'.format(speed_init_type, m, M))
        
        charge = available_particles[type_]['charge']
        mass = available_particles[type_]['mass']
        effective_diameter = available_particles[type_]['effective diameter']

        # parameters of the maxwellian distribution
        # https://stackoverflow.com/questions/63300833/maxwellian-distribution-in-python-scipy
        if(speed_init_type == 'gaussian'): 
            sigma = M
            mu = m
        elif(speed_init_type == 'maxwellian'):
            sigma = M
            mu = m
            a = sigma * np.sqrt(np.pi/(3.0*np.pi - 8.0))
            m_ = 2.0*a*np.sqrt(2.0/np.pi)
            loc = mu - m_
        elif(speed_init_type == 'uniform' or speed_init_type == 'uniform_norm'):
            # uniform distribution parameters
            min_speed_uniform_distribution = m
            max_speed_uniform_distribution = M
        else :
            print("/!\\ Unknown speed init type /!\\")
            return
        mean = 0
        if(debug):print('Number : {}'.format(number_))
        for k in range(number_):
            my_speed = 0

            if(speed_init_type=='gaussian'):
                vx = norm.rvs(mu, sigma)
                vy = norm.rvs(mu, sigma)
                vz = norm.rvs(mu, sigma)
                
            elif(speed_init_type == 'maxwellian'):
                norm_speed = float(maxwell.rvs(loc, a))
                
                theta = random()*2*np.pi
                cTheta = float(np.cos(theta))
                sTheta = float(np.sin(theta))

                phi = random()*2*np.pi
                cPhi = float(np.cos(phi))
                sPhi = float(np.sin(phi))

                vx, vy, vz = norm_speed*sTheta*cPhi, norm_speed*sTheta*sPhi, norm_speed*cTheta

            elif(speed_init_type=='uniform'):
                # norm_speed = min_speed_uniform_distribution+random()*\
                    # (max_speed_uniform_distribution-min_speed_uniform_distribution)
                # direction of the speed
                vx = min_speed_uniform_distribution+random()*(max_speed_uniform_distribution-min_speed_uniform_distribution)
                vy = min_speed_uniform_distribution+random()*(max_speed_uniform_distribution-min_speed_uniform_distribution)
                vz = min_speed_uniform_distribution+random()*(max_speed_uniform_distribution-min_speed_uniform_distribution)

            elif(speed_init_type=='uniform_norm'):
                norm_speed = min_speed_uniform_distribution+random()*\
                                    (max_speed_uniform_distribution-min_speed_uniform_distribution)
                theta = random()*2*np.pi
                cTheta = float(np.cos(theta))
                sTheta = float(np.sin(theta))

                phi = random()*2*np.pi
                cPhi = float(np.cos(phi))
                sPhi = float(np.sin(phi))

                vx, vy, vz = norm_speed*sTheta*cPhi, norm_speed*sTheta*sPhi, norm_speed*cTheta
                
            # my_speed = MyVector(norm_speed*cTheta,norm_speed*sTheta,0)
            my_speed = MyVector(vx, vy, vz)

            x, y = get_correct_initial_positions(zone, offsets, space_size)
            list_particles.append(Particule(charge = charge, radius = effective_diameter/2.0, 
                    mass = mass, part_type = type_, \
                        speed=my_speed, \
                            pos=MyVector(x,y,0), \
                                verbose = verbose))
            mean += my_speed.norm()
            if(debug): 
                print(my_speed)

        if(debug): print("Mean speed init : {} m/s".format(round(mean/number_,2)))

    return np.array(list_particles), mean/number_ # we don't shuffle right now
# TODO : make it work for several particles types ...

def init_particles_flux(wall, direction, nb_particles_to_inject, particles_types, temperatures, drifts):
    list_particles=[]
    e = 1.6e-19 # C
    x1, y1, x2, y2 = wall
    lenght = np.sqrt((x2-x1)**2+(y2-y1)**2)
    for k in range(len(nb_particles_to_inject)):
        type_ = particles_types[k]
        number_ = nb_particles_to_inject[k]
        temperature, drift = temperatures[k], drifts[k]
        
        charge = available_particles[type_]['charge']
        mass = available_particles[type_]['mass']
        effective_diameter = available_particles[type_]['effective diameter']

        sigma = get_gaussian_params_maxwellian(temperature, available_particles[type_]['mass'])
        mu = 0
        for k in range(number_):
            my_speed = 0

            v_drift = drift*direction
            #     
            # vx = norm.rvs(mu, sigma)+v_drift.x # in theory there should be a drift
            # vy = norm.rvs(mu, sigma)+v_drift.y
            # while(direction.inner(MyVector(vx,vy))<=0):
            #     vx = norm.rvs(mu, sigma)+v_drift.x
            #     vy = norm.rvs(mu, sigma)+v_drift.y
            vx, vy = 0, 0
            if(v_drift.y == 0.0):
                vx = sigma * np.sqrt(-2*np.log10((1-random()))) + v_drift.x
                vy = norm.rvs(mu, sigma)
            if(v_drift.x == 0.0):
                vx = norm.rvs(mu, sigma)
                vy = sigma * np.sqrt(-2*np.log10((1-random()))) + v_drift.y
            vz = norm.rvs(mu, sigma)
            my_speed = MyVector(vx, vy, vz)

            # TODO: make something better here
            x, y = x1+random()*(x2-x1)+direction.x*0.01*lenght, y1+random()*(y2-y1)+direction.y*0.01*lenght
            list_particles.append(Particule(charge = charge, radius = effective_diameter/2.0, 
                    mass = mass, part_type = type_, \
                        speed=my_speed, \
                            pos=MyVector(x,y,0)))
    #pprint(" - ".join([part.to_string() for part in list_particles]))
    return np.array(list_particles)

def get_correct_initial_positions(zone, offsets, space_size):
    # same trouble as creating a grid for the space
    # how do we do that best
    offset_x, offset_y = offsets[0], offsets[1]
    lx, ly = space_size[0],space_size[1]
    
    x=random()*lx - offset_x
    y=random()*ly - offset_y
    if(zone!=None):
        while(not zone.inside(Point(x,y))):
            x=random()*lx - offset_x
            y=random()*ly - offset_y
    return x,y
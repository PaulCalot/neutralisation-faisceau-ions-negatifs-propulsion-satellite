from pprint import pprint
import numpy as np
from dolfin import Point
from tqdm import tqdm
from shutil import copyfile, copy 
from os import mkdir
import matplotlib.pyplot as plt
# importing functions and modules useful for launching simulations
from main import get_options, get_min_mean_free_path, init_particles_in_system
from main import square, thruster
from main import CollisionHandler, DSMC
from main import scipy_integrate_solve_ivp, rk4, euler_explicit
from main import DataAnalyser, DataSaver
from main import available_particles
from main import post_processing
from main import convert_list_to_string
from collections import OrderedDict

def simulate(system_cfg_path, simulation_cfg_path, processing_cfg_path, verbose = False):
    system_cfg = system_cfg_path
    simulation_cfg = simulation_cfg_path
    processing_cfg = processing_cfg_path
    
    list_cfg_files = [system_cfg, simulation_cfg, processing_cfg]
    cfg_path_dict = {
        'system':system_cfg,
        'simulation':simulation_cfg,
        'processing':processing_cfg
    }

    options = get_options(cfg_path_dict)
    if(verbose):
        pprint(options)
    system, zone, offset = init_system(options['system'])
    f, args = get_update_fn(system, options['system']['system_type'])

    # TODO there could be an initialization of particles of the system, in which case mean_speed would not be none
    mean_speed=None

    simu = options['simulation']
    integration_scheme = get_scheme(simu)

    system.init_dsmc(mean_speed = mean_speed, integration_scheme = integration_scheme, f = f, f_args = args)
    # simu params
    params_dict, data_analyser = get_saving_dict(options, system)
    
    if(verbose):
        pprint(params_dict)
        system.plot()
        plt.show()
    # loop
    handler = system.dsmc
    MAX_INTEGRATION_STEP =  simu['number_of_steps']
    saving_period = options['processing']['period']
    saving_offset = options['processing']['saving_offset']
    save_test = options['processing']['save']
    dir_saving = options['processing']['path']
    dt = simu['dt']*float(params_dict['min_mean_free_time'])
    t = 0
    for k in tqdm(range(MAX_INTEGRATION_STEP)): #
        system.step(dt, t)
        t+=dt
        list_particles = system.get_list_particles()
        if(verbose and k%saving_period==0): 
            print('Iteration {} : {} particles.'.format(k, len(list_particles)))
            system.plot()
        if(save_test and k > saving_offset and k%saving_period==0 or (k == MAX_INTEGRATION_STEP-1 and k%saving_period!=0)):
            if(len(list_particles) > 0):
                data_analyser.save_everything_to_one_csv(list_particles = list_particles, iteration = k, erase = True)
            
            handler.save_collisions_matrix(name = dir_saving/("test_"+params_dict['id_test']+"_collision_matrix.txt"), iteration = k)
            
    # last saving
    if(save_test):
        update_params(handler, params_dict, data_analyser, dir_saving)
        post_processing(options['processing'])
        for file in list_cfg_files:
            copy(src = file, dst = dir_saving/"cfg_files/")

def processing_only(system_cfg_path, simulation_cfg_path, processing_cfg_path, verbiose = False):
    system_cfg = system_cfg_path
    simulation_cfg = simulation_cfg_path
    processing_cfg = processing_cfg_path

    list_cfg_files = [system_cfg, simulation_cfg, processing_cfg]

    cfg_path_dict = {
        'system':system_cfg,
        'simulation':simulation_cfg,
        'processing':processing_cfg
    }
    options = get_options(cfg_path_dict)
    dir_saving = options['processing']['path']
    
    post_processing(options['processing'])
    for file in list_cfg_files:
        copy(src = file, dst = dir_saving/"cfg_files/")

def init_system(system_options):
    system_type = system_options['system_type']

    system, zone, offset = None, None, None

    if(system_type=='square'):
        system = square(system_options)
        zone = None
        offset = [0,0]
        
    elif(system_type=='thruster'):
        system = thruster(system_options)
        zone = system.get_zone()
        offset = system.get_offset()        
    return system, zone, offset

def get_update_fn(system, system_type):

    if(system_type=='thruster'):
        def f(Y,t,m,q):
            vx=Y[3]
            vy=Y[4]
            vz=Y[5]

            ax = 0 
            ay = 0
            az = 0
            return np.array([vx, vy, vz, ax, ay, az])

        args = []
        # def f(Y,t,m,q,zone,E):
        #     '''
        #     Renvoie la dérivée de Y en y (en faisant un bilan des forces notamment) pr être entré dans RK4
        #     Y=[x, y, z, vx, vy, vz] ce n'est pas le Y de liste_Y
        #     '''
        #     Ex, Ey = E.split(deepcopy=True)
        #     vx=Y[3]
        #     vy=Y[4]
        #     vz=Y[5]
        #     if zone.inside(Point(Y[0],Y[1])):
        #         ax = (1/m) * q * Ex(Y[0], Y[1])
        #         ay = (1/m) * q * Ey(Y[0], Y[1])
        #     else :
        #         ax = 0  #utile si les ki sont hors du mesh,
        #         ay = 0
        #     az=0
        #     return np.array([vx, vy, vz, ax, ay, az])
        
        # args = [system.get_zone(), system.get_E()]
        
    elif(system_type=='square'):
        def f(Y,t,m,q):
            vx=Y[3]
            vy=Y[4]
            vz=Y[5]

            ax = 0 
            ay = 0
            az = 0
            return np.array([vx, vy, vz, ax, ay, az])
        
        args = []
    return f, args


def get_scheme(simu_options):
    integration_scheme = None
    if(simu_options['scheme']=='euler_explicit'):
        integration_scheme = euler_explicit # scipy_integrate_solve_ivp , rk4 , euler_explicit
        # TODO : make a function get_scheme(str)
    return integration_scheme

def get_saving_dict(options, system):
    types = options['system']['particles_types']
    charges = [available_particles[type_]['charge'] for type_ in types]
    masses = [available_particles[type_]['mass'] for type_ in types]
    effective_diameters = [available_particles[type_]['effective diameter'] for type_ in types]

    min_mean_free_path = get_min_mean_free_path(effective_diameters = effective_diameters, particles_densities = options['system']['particles_densities'])
    v_mean = system.get_mean_speed()
    min_mean_free_time = min_mean_free_path/v_mean
    dt = options['simulation']['dt']*min_mean_free_time
    
    MAX_INTEGRATION_STEP = options['simulation']['number_of_steps']

    save_test = options['processing']['save']
    lx, ly, lz = system.get_size()
    offset = system.get_offset()
    res_x, res_y = system.get_resolutions()
    
    saving_period = options['processing']['period']
    dir_saving = options['processing']['path']
    path_saving = dir_saving/(options['processing']['id_test']+'.csv')

    # DSMC params
    N_particles_real = [int(density*system.get_volume()) for density in options['system']['particles_densities']] # this is the REAL number of particles
    N_particles_simu = [int(nb*system.get_number_of_cells()) for nb in options['system']['particles_mean_number_per_cell']]
    Ne = system.get_Ne()

    # dict
    params_dict = OrderedDict()
    
    # general :
    params_dict['id_test']= options['processing']['id_test']
    params_dict['path_to_data'] = str(path_saving) # absolute path

    # on results of the simulation
    params_dict['number_of_collisions'] = '0'
    params_dict['mean_acceptance_rate'] = '0'
    params_dict['mean_proba'] = '0'
    params_dict['mean_vr_norm'] = '0'
    # on system =
    params_dict['nb_cells'] = str(system.get_number_of_cells())
    params_dict['volume'] = str(system.get_volume())
    params_dict['x_min'] = str(offset[0])
    params_dict['x_max'] = str(offset[0]+lx)
    params_dict['y_min'] = str(offset[1])
    params_dict['y_max'] = str(offset[1]+ly)
    # on particles =
    params_dict['total_number_of_particles'] = str(system.get_number_of_particles())
    params_dict['particles_types'] = convert_list_to_string(types)
    params_dict['particles_charges'] = convert_list_to_string(charges)
    params_dict['particles_masses'] = convert_list_to_string(masses)
    params_dict['particles_effective_diameters'] = convert_list_to_string(effective_diameters)
    params_dict['particles_densities'] = convert_list_to_string(options['system']['particles_densities']),
    params_dict['nb_parts_per_type'] = convert_list_to_string(N_particles_simu)
    params_dict['nb_parts_real_per_type'] = convert_list_to_string(N_particles_real)
    params_dict['Ne_per_type'] = convert_list_to_string(Ne)
    params_dict['mean_particles_number_per_cell_per_type'] = convert_list_to_string(options['system']['particles_mean_number_per_cell'])
        # deduction
    params_dict['min_mean_free_path'] = str(min_mean_free_path)
    params_dict['min_mean_free_time'] = str(min_mean_free_time)
        # on speed
    params_dict['speed_init_type'] = convert_list_to_string(options['system']['speed_type'])
    params_dict['v_mean'] = str(v_mean)
    params_dict['vr_max_init'] = str(2*v_mean)
    params_dict['vr_max_final'] = '0'

    # on simulation params =
    params_dict['dt'] = str(dt)
    params_dict['MAX_INTEGRATION_STEP'] = str(MAX_INTEGRATION_STEP)
    params_dict['integration_scheme'] = system.get_integration_scheme().__name__
    #'eta'= eta
    #'loss_charge_proba' = p
        
    # on saving params
    params_dict['saving_period'] = str(saving_period) 

    data_analyser = None
    if(save_test):    
        data_analyser = DataSaver(name_test = options['processing']['id_test']+'.csv', saving_directory = dir_saving)
        data_analyser.save_test_params(dir_saving/'params.csv', params_dict, use_saving_directory = True)
    return params_dict, data_analyser

def update_params(handler, params_dict, data_analyser, path_dir):
    number_of_collisions = handler.get_collisions_count()
    mean_acceptance_rate = handler.get_acceptance_rate()
    mean_vr_norm = handler.get_mean_vr_norm()
    vr_max_final = handler.get_vr_norm()
    mean_proba = handler.get_mean_proba()

    # saving again to the csv.
    params_dict['mean_proba']=str(mean_proba)
    params_dict['number_of_collisions']=str(number_of_collisions)
    params_dict['mean_acceptance_rate']=str(mean_acceptance_rate)
    params_dict['mean_vr_norm']=str(mean_vr_norm)
    params_dict['vr_max']=str(vr_max_final)
    data_analyser.update_saved_params(path_dir/'params.csv', params_dict, use_saving_directory = False)


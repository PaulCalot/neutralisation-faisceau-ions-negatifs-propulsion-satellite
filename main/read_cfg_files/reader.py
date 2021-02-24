import configparser as ConfigParser
#argparse doc: https://docs.python.org/2/library/argparse.html#module-argparse
# https://docs.python.org/3/library/configparser.html#ConfigParser.RawConfigParser.optionxform
import argparse
from main import cfg_tools
# 
from pathlib import Path, PurePath
from os import makedirs
from os.path import isfile, exists

# ----------- get_options ----------------- #

def get_options(cfg_file_paths, verbose = True):
    
    system_options = read_conf_system(cfg_file_paths['system'])
    simulation_options = read_conf_simulation(cfg_file_paths['simulation'])
    processing_options = read_conf_processing(cfg_file_paths['processing'])
    
    # system options
    system_type = str(system_options.system_type)
    speed_type = list(map(str,system_options.speed_type.split(',')))
    speed_param1 = list(map(float,system_options.speed_param1.split(',')))
    speed_param2 = list(map(float,system_options.speed_param2.split(',')))
    particles_types = list(map(cfg_tools.read_args_multiple_types,system_options.particles_types.split(',')))
    particles_mean_number_per_cell = list(map(int,system_options.particles_mean_number_per_cell.split(',')))
    particles_densities = list(map(float,system_options.particles_densities.split(',')))
    #particles_radius = list(map(float,main_options.particles_radius.split(',')))

    flux_in_wall = cfg_tools.read_args_multiple_types(system_options.flux_in_wall)
    flux_out_wall = cfg_tools.read_args_multiple_types(system_options.flux_out_wall)
    temperatures = list(map(cfg_tools.read_args_multiple_types,system_options.temperatures.split(',')))
    drifts = list(map(cfg_tools.read_args_multiple_types,system_options.drifts.split(',')))

    system = {
        'system_type': system_type,
        'speed_type': speed_type,
        'speed_param1': speed_param1,
        'speed_param2': speed_param2,
        'particles_types': particles_types,
        'particles_mean_number_per_cell': particles_mean_number_per_cell,
        'particles_densities': particles_densities,
        #'particles_radius': particles_radius,
        'flux_in_wall':flux_in_wall,
        'flux_out_wall':flux_out_wall,
        'temperatures':temperatures,
        'drifts':drifts

    }

    # system options
    if(system_type == 'square'):
        other_ = get_options_square(system_options)
    elif(system_type == 'thruster'):
        other_ = get_options_thruster(system_options)
    else : 
        if verbose : print("{} system type not recognized. Stopping execution now.".format(system_type))
        exit()
    system.update(other_)

    # processing options
        # saving options
    save = cfg_tools.str_to_bool(processing_options.save)
    period = int(processing_options.period)
    saving_offset = int(processing_options.saving_offset)
    path = Path(__file__).parent.parent.parent.absolute() / \
         'results' / str(processing_options.path)
    id_test = str(processing_options.id_test) # list(map(str, processing_options.ids_test.split(',')))

    if(save and not exists(path)):
        makedirs(path)
        if(not exists(path/'figures/')) : makedirs(path/'figures/')
        if(not exists(path/'cfg_files/')) : makedirs(path/'cfg_files/')
    
        # post_processing options
    compute_system_evolution = cfg_tools.str_to_bool(processing_options.compute_system_evolution)
    compute_hist_distribution_evolution = cfg_tools.str_to_bool(processing_options.compute_hist_distribution_evolution)
    compute_spatial_distribution = cfg_tools.str_to_bool(processing_options.compute_spatial_distribution)
    compute_temperature = cfg_tools.str_to_bool(processing_options.compute_temperature)
    frames_to_compute = list(map(cfg_tools.read_args_multiple_types, \
         processing_options.frames_to_compute.split(',')))
    compute_collisions = cfg_tools.str_to_bool(processing_options.compute_collisions)
    compute_density = cfg_tools.str_to_bool(processing_options.compute_density)
    args_density = list(map(int, \
         processing_options.args_density.split(',')))
    compute_number_of_particles = cfg_tools.str_to_bool(processing_options.compute_number_of_particles)

    
    merge_csv = cfg_tools.str_to_bool(processing_options.merge_csv)
    files_to_merge = list(map(str, processing_options.files_to_merge.split(',')))

    processing = {
        'save' : save,
        'period':period,
        'saving_offset':saving_offset,
        'path':path,
        'id_test':id_test,
        'compute_system_evolution' : compute_system_evolution,
        'compute_hist_distribution_evolution' :compute_hist_distribution_evolution,
        'compute_spatial_distribution':compute_spatial_distribution,
        'compute_temperature':compute_temperature,
        'compute_density':compute_density,
        'args_density':args_density,
        'compute_number_of_particles':compute_number_of_particles,
        'frames_to_compute':frames_to_compute,
        'compute_collisions':compute_collisions,
        'merge_csv':merge_csv,
        'files_to_merge':files_to_merge,
    }

    # simulation options
    seed = cfg_tools.read_args_multiple_types(simulation_options.seed) # can also be None
    scheme = str(simulation_options.scheme)
    dt = float(simulation_options.dt)
    number_of_steps = int(simulation_options.number_of_steps)

    simulation = {
        'seed':seed,
        'scheme':scheme,
        'dt':dt,
        'number_of_steps':number_of_steps
    }

    # final options
    options = {
        'system':system,
        'simulation':simulation,
        'processing':processing
    }
    
    return options

def get_options_square(system_options):

    resolution = list(map(int,system_options.resolution.split(',')))
    size = list(map(float,system_options.size.split(',')))
    lz = float(system_options.lz) 

    return {
        'resolution':resolution,
        'size':size,
        'lz':lz
    }
    
def get_options_thruster(system_options):

    resolution = list(map(int,system_options.resolution.split(',')))
    lz = float(system_options.lz) 

    general = {
        'resolution':resolution,
        'lz':lz
    }

    L_mot = float(system_options.L_mot)
    l_mot = float(system_options.l_mot)
    L_1 = float(system_options.L_1)
    l_1 = float(system_options.l_1)
    L_2 = float(system_options.L_2)
    l_2 = float(system_options.l_2)
    delta_vert_12 = float(system_options.delta_vert_12)
    L_vacuum = float(system_options.L_vacuum)
    l_vacuum = float(system_options.l_vacuum)
    mesh_resolution = int(system_options.mesh_resolution)
    refine_mesh = bool(system_options.refine_mesh)

    geometry = {
        'L_mot':L_mot,
        'l_mot':l_mot,
        'L_1':L_1,
        'l_1':l_1,
        'L_2':L_2,
        'l_2':l_2,
        'delta_vert_12':delta_vert_12,
        'L_vacuum':L_vacuum,
        'l_vacuum':l_vacuum,
        'mesh_resolution':mesh_resolution,
        'refine_mesh':refine_mesh,
    }

    phi_top_mot = cfg_tools.read_args_multiple_types(system_options.phi_top_mot)
    phi_bord_mot = cfg_tools.read_args_multiple_types(system_options.phi_bord_mot)
    phi_electrode1 = cfg_tools.read_args_multiple_types(system_options.phi_electrode1)
    phi_inter_electrode = cfg_tools.read_args_multiple_types(system_options.phi_inter_electrode)
    phi_electrode2 = cfg_tools.read_args_multiple_types(system_options.phi_electrode2)
    phi_sup_vacuum = cfg_tools.read_args_multiple_types(system_options.phi_sup_vacuum)
    phi_inf_vacuum = cfg_tools.read_args_multiple_types(system_options.phi_inf_vacuum)
    rho_elec = float(system_options.rho_elec) 

    phi = {
        'phi_top_mot':phi_top_mot,
        'phi_bord_mot' :phi_bord_mot,
        'phi_electrode1' :phi_electrode1,
        'phi_inter_electrode':phi_inter_electrode,
        'phi_electrode2' :phi_electrode2,
        'phi_sup_vacuum' :phi_sup_vacuum,
        'phi_inf_vacuum':phi_inf_vacuum,
        'rho_elec':rho_elec
    }
    
    return {
        'general':general,
        'phi':phi,
        'geometry':geometry
    }
# ------------------ Reader cfg files ------------------ #
# TODO : add the other systems (only 'square' for now)

def read_conf_system(cfg_file_path):
    # Initializing dummy class with cfg folder path
    options = cfg_tools.Options(cfg_file_path)

    # Reading the config file with config parser
    config = ConfigParser.RawConfigParser()
    config.optionxform = lambda option: option
    #config = ConfigParser.ConfigParser()
    config.read(options.cfg)

    if(config.has_section('square')):
        read_conf_square(options, config)
        options.system_type = 'square'
    elif(config.has_section('thruster')):
        read_conf_thruster(options, config)
        options.system_type = 'thruster'
    else:
        print('System not present')
        return

    #[speed]
    options.speed_type = config.get('speed','speed_type')
    options.speed_param1 = config.get('speed','speed_param1')
    options.speed_param2 = config.get('speed','speed_param2')

    #[particles]
    options.particles_types = config.get('particles','particles_types')
    options.particles_mean_number_per_cell = config.get('particles','particles_mean_number_per_cell')
    options.particles_densities = config.get('particles','particles_densities')

    #[collisions]
    #options.eta = config.get('collisions','eta')
    #options.rho = config.get('collisions','rho')

    #[flux]
    if(config.has_section('flux')):
        options.flux_in_wall = config.get('flux','flux_in_wall')
        options.flux_out_wall = config.get('flux','flux_out_wall')
        options.temperatures = config.get('flux','temperatures')
        options.drifts = config.get('flux','drifts')
    else:
        options.flux_in_wall = 'None'
        options.flux_out_wall = 'None'
        options.temperatures = 'None'
        options.drifts = 'None'
    return options

def read_conf_simulation(cfg_file_path):
    # Initializing dummy class with cfg folder path
    options = cfg_tools.Options(cfg_file_path)

    # Reading the config file with config parser
    config = ConfigParser.ConfigParser()
    config.read(options.cfg)

    options.seed = config.get('params','seed')
    options.scheme = config.get('params','scheme')
    options.dt = config.get('params','dt')
    options.number_of_steps = config.get('params','number_of_steps')

    return options

def read_conf_processing(cfg_file_path):
    options = cfg_tools.Options(cfg_file_path)

    # Reading the config file with config parser
    config = ConfigParser.ConfigParser()
    config.read(options.cfg)

    options.save = config.get('general','save')
    options.postprocess = config.get('general','postprocess')
    options.period = config.get('general','period')
    if(config.has_option('general', 'saving_offset')):
        options.saving_offset = config.get('general','saving_offset')
    else:
        options.saving_offset = 0
    options.path = config.get('general','path')
    options.id_test = config.get('general','id_test')

    options.compute_system_evolution = config.get('processing options','compute_system_evolution')
    options.compute_hist_distribution_evolution = config.get('processing options','compute_hist_distribution_evolution')
    options.compute_spatial_distribution = config.get('processing options','compute_spatial_distribution')
    options.compute_temperature = config.get('processing options','compute_temperature')
    options.frames_to_compute = config.get('processing options','frames_to_compute')
    options.compute_collisions = config.get('processing options','compute_collisions')
    options.compute_density = config.get('processing options','compute_density')
    options.args_density = config.get('processing options','args_density')
    options.compute_number_of_particles = config.get('processing options','compute_number_of_particles')

    options.merge_csv = config.get('merge file','merge_csv')
    options.files_to_merge = config.get('merge file','files_to_merge')

    return options

def read_conf_square(options, config):
    # Reading the config file with config parser

    options.resolution = config.get('square','resolution')
    options.size = config.get('square','size')
    options.lz = config.get('square','lz')

    return options

def read_conf_thruster(options, config):

    options.resolution = config.get('thruster','resolution')
    options.lz = config.get('thruster','lz')

    options.L_mot = config.get('thruster','L_mot')
    options.l_mot = config.get('thruster','l_mot')
    options.L_1 = config.get('thruster','L_1')
    options.l_1 = config.get('thruster','l_1')
    options.L_2 = config.get('thruster','L_2')
    options.l_2 = config.get('thruster','l_2')
    options.delta_vert_12 = config.get('thruster','delta_vert_12')
    options.L_vacuum = config.get('thruster','L_vacuum')
    options.l_vacuum = config.get('thruster','l_vacuum')
    options.mesh_resolution = config.get('thruster','mesh_resolution')
    options.refine_mesh = config.get('thruster','refine_mesh')

    options.phi_top_mot = config.get('thruster','phi_top_mot')
    options.phi_bord_mot = config.get('thruster','phi_bord_mot')
    options.phi_electrode1 = config.get('thruster','phi_electrode1')
    options.phi_inter_electrode = config.get('thruster','phi_inter_electrode')
    options.phi_electrode2 = config.get('thruster','phi_electrode2')
    options.phi_sup_vacuum = config.get('thruster','phi_sup_vacuum')
    options.phi_inf_vacuum = config.get('thruster','phi_inf_vacuum')
    options.rho_elec = config.get('thruster','rho_elec')

    return options

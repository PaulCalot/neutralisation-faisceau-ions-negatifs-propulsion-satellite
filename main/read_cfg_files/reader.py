import configparser as ConfigParser
#argparse doc: https://docs.python.org/2/library/argparse.html#module-argparse

import argparse
from main import cfg_tools
# 
from pathlib import Path, PurePath

# ----------- get_options ----------------- #

def get_options(cfg_file_paths, verbose = True):
    
    main_options = read_conf_main(cfg_file_paths['main'])
    saving_options = read_conf_saving(cfg_file_paths['saving'])
    simulation_options = read_conf_simulation(cfg_file_paths['simulation'])
    
    # main options
    system_type = str(main_options.system_type)
    id_test = int(main_options.id_test)
    speed_type = list(map(str,main_options.speed_type.split(',')))
    speed_param1 = list(map(float,main_options.speed_param1.split(',')))
    speed_param2 = list(map(float,main_options.speed_param2.split(',')))
    particles_types = list(map(str,main_options.particles_types.split(',')))
    particles_mean_number_per_cell = list(map(int,main_options.particles_mean_number_per_cell.split(',')))
    particles_densities = list(map(float,main_options.particles_densities.split(',')))
    particles_radius = list(map(float,main_options.particles_radius.split(',')))

    main = {
        'system_type': system_type,
        'id_test' : id_test,
        'speed_type': speed_type,
        'speed_param1': speed_param1,
        'speed_param2': speed_param2,
        'particles_types': particles_types,
        'particles_mean_number_per_cell': particles_mean_number_per_cell,
        'particles_densities': particles_densities,
        'particles_radius': particles_radius,
    }

    # saving options
    period = int(saving_options.period)
    path = Path(__file__).parent.parent.parent.absolute() / 'results' / str(saving_options.path)
    # absolute path to the saving directory
    saving = {
        'period':period,
        'path':path
    }

    # simulation options
    scheme = str(simulation_options.scheme)
    dt = float(simulation_options.dt)
    number_of_steps = int(simulation_options.number_of_steps)

    simulation = {
        'scheme':scheme,
        'dt':dt,
        'number_of_steps':number_of_steps
    }

    options = {
        'main':main,
        'saving':saving,
        'simulation':simulation
    }
    # system options
    if(system_type == 'square'):
        options['square'] = get_options_square(cfg_file_paths['system'])
    else : 
        if verbose : print("{} system type not recognized. Stopping execution now.".format(system_type))
        exit()
    
    return options

def get_options_square(cfg_file_path):
    system_options = read_conf_square(cfg_file_path)

    factor = float(system_options.factor)
    size = list(map(int,system_options.size.split(',')))
    lz = float(system_options.lz) 

    return {
        'factor':factor,
        'size':size,
        'lz':lz
    }

# ------------------ Reader cfg files ------------------ #
# TODO : add the other systems (only 'square' for now)

def read_conf_main(cfg_file_path):
    # Initializing dummy class with cfg folder path
    options = cfg_tools.Options(cfg_file_path)

    # Reading the config file with config parser
    Config = ConfigParser.ConfigParser()
    Config.read(options.cfg)

    #[id tes]
    options.id_test = Config.get('id_test','id_test')
    
    #[system]
    options.system_type = Config.get('system','system_type')

    #[speed]
    options.speed_type = Config.get('speed','speed_type')
    options.speed_param1 = Config.get('speed','speed_param1')
    options.speed_param2 = Config.get('speed','speed_param2')

    #[particles]
    options.particles_types = Config.get('particles','particles_types')
    options.particles_mean_number_per_cell = Config.get('particles','particles_mean_number_per_cell')
    options.particles_densities = Config.get('particles','particles_densities')
    options.particles_radius = Config.get('particles','particles_radius')

    #[collisions]
    #options.eta = Config.get('collisions','eta')
    #options.rho = Config.get('collisions','rho')

    return options

def read_conf_saving(cfg_file_path):
    # Initializing dummy class with cfg folder path
    options = cfg_tools.Options(cfg_file_path)

    # Reading the config file with config parser
    Config = ConfigParser.ConfigParser()
    Config.read(options.cfg)

    options.period = Config.get('params','period')
    options.path = Config.get('params','path')

    return options

def read_conf_simulation(cfg_file_path):
    # Initializing dummy class with cfg folder path
    options = cfg_tools.Options(cfg_file_path)

    # Reading the config file with config parser
    Config = ConfigParser.ConfigParser()
    Config.read(options.cfg)

    options.scheme = Config.get('params','scheme')
    options.dt = Config.get('params','dt')
    options.number_of_steps = Config.get('params','number_of_steps')

    return options

def read_conf_square(cfg_file_path):
    # Initializing dummy class with cfg folder path
    options = cfg_tools.Options(cfg_file_path)

    # Reading the config file with config parser
    Config = ConfigParser.ConfigParser()
    Config.read(options.cfg)

    options.factor = Config.get('system','factor')
    options.size = Config.get('system','size')
    options.lz = Config.get('system','lz')

    return options


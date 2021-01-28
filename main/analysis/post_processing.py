from main import merge_tests_summary
from main import DataAnalyser
import numpy as np

def post_processing(options):

    path_to_file = options['path']
    ids_test = options['ids_test'] # [1,2,3,4] # [2] # [1,2,3,4] #[1,2,3,4]

    compute_system_evolution = options['compute_system_evolution']
    compute_hist_distribution_evolution = options['compute_hist_distribution_evolution']
    compute_spatial_distribution = options['compute_spatial_distribution']
    compute_temperature = options['compute_temperature']
    frames_to_compute = options['frames_to_compute']
    merge_csv = options['merge_csv']

    files_to_merge = options['file_to_merge']

    if(merge_csv):
        merge_tests_summary(files_to_merge, path_to_file)

    for indx, id_test in enumerate(ids_test):
        data_analyser = DataAnalyser(path_to_file)
        data_analyser.load_test(id_test)
        mass_particles = data_analyser.get('mass') # will return, for the current test, the value of the column 'mass' if it exists
        period = data_analyser.get('dt')*data_analyser.get('MAX_INTEGRATION_STEP')
        xy_min = data_analyser.get('xy_min')
        xy_max = data_analyser.get('xy_max')
        max_speed = data_analyser.get('vr_max')
        mean_particles_number_per_cell = data_analyser.get('mean_particles_number_per_cell')
        nb_cells = data_analyser.get('nb_cells')
        
        if(compute_system_evolution):
            data_analyser.draw_particles()

        if(compute_hist_distribution_evolution):
            # speed_norm , speed_norm_squared , vz
            data_analyser.draw_hist_distribution('vx',plot_gaussian=True)
            data_analyser.draw_hist_distribution('vy',plot_gaussian=True)
            data_analyser.draw_hist_distribution('speed_norm', plot_maxwellian = True)
            data_analyser.draw_hist_distribution('vz', plot_gaussian=True)

        if(compute_spatial_distribution):
            data_analyser.draw_spatial_distribution(None, vmin = 0, vmax = 2*mean_particles_number_per_cell) 
            data_analyser.draw_spatial_distribution('vx', vmin = -max_speed, vmax = max_speed)
            data_analyser.draw_spatial_distribution('vy', vmin = -max_speed, vmax = max_speed)
            data_analyser.draw_spatial_distribution('vz', vmin = -max_speed, vmax = max_speed)
            data_analyser.draw_spatial_distribution('speed_norm', vmin = 0.01*max_speed**2, vmax = max_speed**2)

        for which in frames_to_compute:
            data_analyser.draw_particles_frame(which = which, save_frame=True, vmin = xy_min, vmax = xy_max)
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vx', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'r')
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vy', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'g')
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'speed_norm', save_frame = True, plot_maxwellian = True, plot_gaussian = False, density = True, range = None, color = 'k')
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'speed_norm_squared', save_frame = True, plot_maxwellian = False, plot_gaussian = False, density = True, range = None, color = 'k')
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vz', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'b')
            data_analyser.draw_spatial_distribution_frame(which, None, grid_size = int(np.sqrt(nb_cells)), vmin = 0, vmax = 2*mean_particles_number_per_cell)

        if compute_temperature:
            data_analyser.draw_Temperature_evolution(period, tau_init = 5e-5, particles_mass= mass_particles, begin = 0.00, end = 1.0) # g/mol

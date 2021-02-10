from main import merge_tests_summary
from main import DataAnalyser
import numpy as np
from main import convert_string_to_list

def post_processing(options):

    path_to_file = options['path']
    ids_test = options['id_test'] # [1,2,3,4] # [2] # [1,2,3,4] #[1,2,3,4]

    compute_system_evolution = options['compute_system_evolution']
    compute_hist_distribution_evolution = options['compute_hist_distribution_evolution']
    compute_spatial_distribution = options['compute_spatial_distribution']
    compute_temperature = options['compute_temperature']
    frames_to_compute = options['frames_to_compute']
    merge_csv = options['merge_csv']
    
    compute_collisions = options['compute_collisions']

    files_to_merge = options['files_to_merge']

    if(merge_csv):
        merge_tests_summary(files_to_merge, path_to_file)

    if(type(ids_test)!=list):
       ids_test = [ids_test]

    for indx, id_test in enumerate(ids_test):
        data_analyser = DataAnalyser(path_to_file/'params.csv')
        data_analyser.load_test(id_test)

         # particles
        particles_types = convert_string_to_list(data_analyser.get_param('particles_types'))
        particles_effective_diameters = convert_string_to_list(data_analyser.get_param('particles_effective_diameters'))
        particles_charges = convert_string_to_list(data_analyser.get_param('particles_charges'))
        particles_masses = convert_string_to_list(data_analyser.get_param('particles_masses')) # will return, for the current test, the value of the column 'mass' if it exists

        mean_particles_number_per_cell_per_type = convert_string_to_list(data_analyser.get_param('mean_particles_number_per_cell_per_type'))
        period = data_analyser.get_param('dt')*data_analyser.get_param('MAX_INTEGRATION_STEP')
        total_number_of_particles_per_cell = sum(mean_particles_number_per_cell_per_type)
        # space limit
        x_min = data_analyser.get_param('x_min')
        x_max = data_analyser.get_param('x_max')
        y_min = data_analyser.get_param('y_min')
        y_max = data_analyser.get_param('y_max')
        
        # speed
        max_speed = data_analyser.get_param('vr_max_final')

        # system
        nb_cells = data_analyser.get_param('nb_cells')
        
        if(compute_system_evolution):
            data_analyser.draw_particles()

        if(compute_hist_distribution_evolution):
            # speed_norm , speed_norm_squared , vz
            data_analyser.draw_hist_distribution('vx',plot_gaussian=True)
            data_analyser.draw_hist_distribution('vy',plot_gaussian=True)
            data_analyser.draw_hist_distribution('speed_norm', plot_maxwellian = True)
            data_analyser.draw_hist_distribution('vz', plot_gaussian=True)

        if(compute_spatial_distribution):
            data_analyser.draw_spatial_distribution(None, vmin = 0, vmax = 2*total_number_of_particles_per_cell.sum()) 
            data_analyser.draw_spatial_distribution('vx', vmin = -max_speed, vmax = max_speed)
            data_analyser.draw_spatial_distribution('vy', vmin = -max_speed, vmax = max_speed)
            data_analyser.draw_spatial_distribution('vz', vmin = -max_speed, vmax = max_speed)
            data_analyser.draw_spatial_distribution('speed_norm', vmin = 0.01*max_speed**2, vmax = max_speed**2)

        for which in frames_to_compute:
            data_analyser.draw_particles_frame(which = which, save_frame=True,  x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max)
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vx', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'r')
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vy', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'g')
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'speed_norm', save_frame = True, plot_maxwellian = True, plot_gaussian = False, density = True, range = None, color = 'k')
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'speed_norm_squared', save_frame = True, plot_maxwellian = False, plot_gaussian = False, density = True, range = None, color = 'k')
            data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vz', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'b')
            #data_analyser.draw_spatial_distribution_frame(which, None, grid_size = int(np.sqrt(nb_cells)), vmin = 0, vmax = 2*total_number_of_particles_per_cell)

        if compute_temperature:
            for mass in particles_masses:
                data_analyser.draw_Temperature_evolution(period, tau_init = 1e-3, particles_mass= mass, begin = 0.00, end = 1.0) # g/mol

        if compute_collisions:
            from main import collision_frequency_th_V, collision_frequency_th_T
            # TODO : make multi particles types
            nb_parts = data_analyser.get_param('x_min')
            
            dt = data_analyser.get_param('dt')
            number_of_dt = data_analyser.get_param('MAX_INTEGRATION_STEP')
            
            Ne = convert_string_to_list(data_analyser.get_param('Ne_per_type'))
            nb_coll = data_analyser.get_param('number_of_collisions')
            number_of_collisions = Ne[0]*nb_coll
            volume = data_analyser.get_param('volume')
            mean_free_path = data_analyser.get_param('min_mean_free_path')

            v_mean = data_analyser.get_param('v_mean')
            particle_density = data_analyser.get_param('particles_densities')

            f_th_ = collision_frequency_th_V(v_mean, mean_free_path, n = particle_density)
            expected_number_of_collision = f_th_*volume*number_of_dt*dt

            print('For test {}, N_e = {:e} vs N_t = {:e}.'.format(id_test , number_of_collisions, expected_number_of_collision))
        
        nb_frames = 298
        bins_x, bins_y = 40,60

        data_analyser.plot_mean_density_evolution_1D(direction = 'y', bins = bins_y, frames = nb_frames) # last 10 frames, with period 10 means we are going to take the last 100 time step... (in theory)
        data_analyser.plot_mean_density_evolution_1D(direction = 'x', bins = bins_x, frames = nb_frames) # last 10 frames, with period 10 means we are going to take the last 100 time step... (in theory)
        data_analyser.plot_mean_density_evolution_by('y', bins = bins_y//2, bins_ = bins_x//2, frames = nb_frames)
        data_analyser.plot_mean_density_evolution_by('x', bins = bins_x//2, bins_ = bins_y//2, frames = nb_frames)

        Ne = convert_string_to_list(data_analyser.get_param('Ne_per_type'))
        Ne_mean = np.mean(np.array(Ne))
        volume = data_analyser.get_param('volume')/(bins_x*bins_y)
        def f(raw):
            return 2.0*Ne_mean / (volume*nb_frames) # I added 2.0 here because somehow plt.hexbin multiply by two the number of boxes...

        data_analyser.draw_spatial_distribution_frame(which = "last", name = ('density',f, np.sum), save_frame = True, grid_size = (bins_x,bins_y), mean_over_frames = nb_frames) 
        data_analyser.draw_spatial_distribution_frame(which = "last", name = 'speed_norm', save_frame = True, grid_size = (bins_x,bins_y), mean_over_frames = nb_frames) 
        data_analyser.draw_spatial_distribution_frame(which = "last", name = 'vx', save_frame = True, grid_size = (bins_x,bins_y), mean_over_frames = nb_frames) 
        data_analyser.draw_spatial_distribution_frame(which = "last", name = 'vy', save_frame = True, grid_size = (bins_x,bins_y), mean_over_frames = nb_frames) 
        data_analyser.draw_spatial_distribution_frame(which = "last", name = 'vz', save_frame = True, grid_size = (bins_x,bins_y), mean_over_frames = nb_frames) 

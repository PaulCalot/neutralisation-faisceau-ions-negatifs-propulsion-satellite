from modules.data_analysis import DataAnalyser, merge_tests_summary

# mass particles
from modules.utils import get_mass_part

merge_csv = False
compute_system_evolution = False
compute_hist_distribution_evolution = False
compute_spatial_distribution = False
compute_temperature = True
frames_to_compute = ["first","last"] # ["first","last"]

ids_test = [1] #[1,2,3] #[1,2,3,4]
periods = [1.17228686999623e-05, 1.17228686999623e-05, 1.17228686999623e-05] # 5.86143434998116e-06, 1.17228686999623e-05]
tests_summary_file_name='gaussian'

if(merge_csv):
    names = ['test_rapport_inter_serie3_1', 'test_rapport_inter_serie3_2','test_rapport_inter_serie3_3','test_rapport_inter_serie3_4']
    output_name = 'test_rapport_inter_serie3'
    merge_tests_summary(names, output_name)

for indx, id_test in enumerate(ids_test):
    data_analyser = DataAnalyser(tests_summary_file_name)
    data_analyser.load_test(id_test)

    if(compute_system_evolution):
        data_analyser.draw_particles()

    if(compute_hist_distribution_evolution):
        # speed_norm , speed_norm_squared , vz
        data_analyser.draw_hist_distribution('vx',plot_gaussian=True)
        data_analyser.draw_hist_distribution('vy',plot_gaussian=True)
        data_analyser.draw_hist_distribution('speed_norm', plot_maxwellian = True)
        data_analyser.draw_hist_distribution('vz', plot_gaussian=True)
    if(compute_spatial_distribution):
        data_analyser.draw_spatial_distribution(None, vmin = 0, vmax = 150)
        data_analyser.draw_spatial_distribution('vx', vmin = -5e3, vmax = 5e3)
        data_analyser.draw_spatial_distribution('vy', vmin = -5e3, vmax = 5e3)
        data_analyser.draw_spatial_distribution('vz', vmin = -5e3, vmax = 5e3)
        data_analyser.draw_spatial_distribution('speed_norm', vmin = 1e6, vmax = 25e6)

    for which in frames_to_compute:
        data_analyser.draw_particles_frame(which = which, save_frame=True, vmin = 0, vmax = 0.042)
        data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vx', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'r')
        data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vy', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'g')
        data_analyser.draw_hist_distribution_frame(which = which, value_name = 'speed_norm', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'k')
        data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vz', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None, color = 'b')
        data_analyser.draw_spatial_distribution_frame(which, None, grid_size = 9, vmin = 0, vmax = 60 )

    if compute_temperature:
        IODINE_MASS = get_mass_part(53, 53, 88) 
        data_analyser.draw_Temperature_evolution(periods[indx], tau_init = 5e-5, particles_mass= IODINE_MASS, begin = 0.00, end = 1.0) # g/mol
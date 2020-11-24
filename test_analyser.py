from modules.data_analysis import DataAnalyser, merge_tests_summary

merge_csv = False
compute_system_evolution = False
compute_hist_distribution_evolution = False
compute_spatial_distribution = False
frames_to_compute = ["first","last"]

id_test = 14
tests_summary_file_name='tests_summary_v4'

if(merge_csv):
    names = ['tests_summary_7', 'tests_summary_8','tests_summary_9',\
        'tests_summary_10','tests_summary_11']
    output_name = 'tests_summary_v2'
    merge_tests_summary(names, output_name)

data_analyser = DataAnalyser(tests_summary_file_name)
data_analyser.load_test(id_test)


if(compute_system_evolution):
    data_analyser.draw_particles()

if(compute_hist_distribution_evolution):
    # speed_norm , speed_norm_squared , vz
    data_analyser.draw_hist_distribution('vx',plot_gaussian=True)
    data_analyser.draw_hist_distribution('vy',plot_gaussian=True)
    data_analyser.draw_hist_distribution('speed_norm_squared', plot_maxwellian = True)
    data_analyser.draw_hist_distribution('vz', plot_gaussian=True)
if(compute_spatial_distribution):
    data_analyser.draw_spatial_distribution(None, vmin = 0, vmax = 150)
    data_analyser.draw_spatial_distribution('vx', vmin = -5e3, vmax = 5e3)
    data_analyser.draw_spatial_distribution('vy', vmin = -5e3, vmax = 5e3)
    data_analyser.draw_spatial_distribution('vz', vmin = -5e3, vmax = 5e3)
    data_analyser.draw_spatial_distribution('speed_norm_squared', vmin = 1e6, vmax = 25e6)

for which in frames_to_compute:
    data_analyser.draw_particles_frame(which = which, save_frame=True)
    data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vx', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None)
    data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vy', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None)
    data_analyser.draw_hist_distribution_frame(which = which, value_name = 'speed_norm_squared', save_frame = True, plot_maxwellian = True, plot_gaussian = True, density = True, range = None)
    data_analyser.draw_hist_distribution_frame(which = which, value_name = 'vz', save_frame = True, plot_maxwellian = False, plot_gaussian = True, density = True, range = None)
    data_analyser.draw_spatial_distribution_frame(which, None, grid_size = 10, vmin = 0, vmax = 60  )

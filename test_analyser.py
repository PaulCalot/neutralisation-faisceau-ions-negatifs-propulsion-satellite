from modules.data_analysis import DataAnalyser


id_test = 4
tests_summary_file_name='tests_summary_2'
data_analyser = DataAnalyser(tests_summary_file_name)
data_analyser.load_test(id_test)

data_analyser.draw_particles()
#speed_norm , speed_norm_squared , vz
#data_analyser.draw_hist_distribution('vx',plot_gaussian=True)
#data_analyser.draw_hist_distribution('vy',plot_gaussian=True)
#data_analyser.draw_hist_distribution('speed_norm_squared')
#data_analyser.draw_hist_distribution('vz', plot_gaussian=True)
#data_analyser.draw_spatial_distribution('vx', vmin = -5e3, vmax = 5e3)
#data_analyser.draw_spatial_distribution('vy', vmin = -5e3, vmax = 5e3)
data_analyser.draw_spatial_distribution('vz', vmin = -5e3, vmax = 5e3)
#data_analyser.draw_spatial_distribution('speed_norm_squared', vmin = 1e6, vmax = 25e6)
#data_analyser.draw_spatial_distribution()

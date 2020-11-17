from modules.data_analysis import DataAnalyser


id_test = 4
tests_summary_file_name='tests_summary_2'
data_analyser = DataAnalyser(tests_summary_file_name)
data_analyser.load_test(id_test)

#data_analyser.draw_particles()
# speed_norm , speed_norm_squared , vz
data_analyser.draw_hist_distribution('vx')
data_analyser.draw_hist_distribution('vy')
data_analyser.draw_hist_distribution('speed_norm_squared')

#data_analyser.draw_spatial_distribution()

from modules.data_analysis import DataAnalyser

id_test = 1
tests_summary_file_name='tests_summary'
data_analyser = DataAnalyser(tests_summary_file_name)
data_analyser.load_test(id_test)
#data_analyser.draw_speed_norm_distribution()
#data_analyser.draw_speed_spatial_distribution()
data_analyser.draw_particles_density_distribution()
#data_analyser.draw_particles()



id_test = 2
tests_summary_file_name='tests_summary'
data_analyser = DataAnalyser(tests_summary_file_name)
data_analyser.load_test(id_test)
#data_analyser.draw_speed_norm_distribution()
#data_analyser.draw_speed_spatial_distribution()
data_analyser.draw_particles_density_distribution()
#data_analyser.draw_particles()


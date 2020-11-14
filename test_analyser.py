from modules.data_analysis import DataAnalyser

id_test = 0
tests_summary_file_name='tests_summary'
data_analyser = DataAnalyser(tests_summary_file_name)
data_analyser.draw_speed_norm_distribution(id_test)
from modules.utils import collision_frequency_th_T, collision_frequency_th_V, collision_frequency_exp
import pandas as pd

csv_file = "tests"
test_to_compute = [1,2,3,4]

df = pd.read_csv(csv_file+'.csv', sep = ',', header=0,index_col=0) # take first line as headers

for i,test in enumerate(test_to_compute):
    test_params = df.loc[test]
    path_to_data = test_params.at['path_to_data'] # df.loc[5].at['B']
    nb_parts = test_params.at['total_number_of_particles']
    
    dt = test_params.at['dt']
    number_of_dt = test_params.at['MAX_INTEGRATION_STEP']

    number_of_collisions = test_params.at['number_of_collisions']*test_params.at['Ne']
    volume = test_params.at['volume']
    mean_free_path = test_params.at['mean_free_path']
    
    f_th_ = collision_frequency_th_V(test_params.at['vmean_exp'],mean_free_path, n = test_params.at['real_particle_density'])
    expected_number_of_collision = f_th_*volume*number_of_dt*dt
    print('For test {}, N_e = {:e} vs N_t = {}.'.format(test , number_of_collisions,expected_number_of_collision))
    #df_data = pd.read_csv(path_to_data+'.csv', sep = ',', header=0, index_col=0) # take first line as headers
    
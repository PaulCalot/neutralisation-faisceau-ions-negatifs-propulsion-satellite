
# what is a callback function : https://stackoverflow.com/questions/824234/what-is-a-callback-function

import csv
from os.path import isfile

# local import
from .particules import Particule

# TODO : add the possibility to do some preprocessing before saving 
# (as it risks to be too much data to save)

class DataSaver :
    def __init__(self, particles, name_test, saving_directory):
        self.particles = particles
        self.name_test = name_test
        self.saving_directory = saving_directory

    def save_to_csv_(self, name, list_data, erase = True, use_saving_directory = True): # by default will erase what was already here
        
        path = self.saving_directory + name if use_saving_directory else name

        if(isfile(path + '.csv') and not erase):
            mode = 'a'
        else :
            mode = 'w'
        with open(path + '.csv', mode = mode) as csv_file:
                csv_writer = csv.writer(csv_file, delimiter=',')
                for data in list_data:
                    csv_writer.writerow(data)
    
    def save_everything_to_multiple_csv(self, iter):
        # iter : iteration number
        list_data = [self.particles[0].get_headers()]
        for k in range(len(self.particles)):
            list_data.append(self.particles[k].to_list())
        self.save_to_csv_(self.name_test + "_iter{}".format(iter), list_data)

    def save_everything_to_one_csv(self):
        # will add to the following
        list_data = [self.particles[0].get_headers()]
        for k in range(len(self.particles)):
            list_data.append(self.particles[k].to_list())
        self.save_to_csv_(self.name_test, list_data, erase = False)

    # saving test params
    # TODO : could make a function that save a dictionnary (could be reusable)
    # TODO : make sure to erase the row (if it exists) that has the same test id (or automatically rename it) (is it possible?)
    def save_test_params(self, tests_summary_file_name, params_dict, use_saving_directory = True):
        path = self.saving_directory + tests_summary_file_name if use_saving_directory else tests_summary_file_name

        if(isfile(path + '.csv')):
            mode = 'a'
        else :
            mode = 'w'
        
        with open(path + '.csv', mode = mode) as csv_file:
                fieldnames = params_dict.keys()
                csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                if(mode == "w"):
                    # first time, we add the keys
                    csv_writer.writeheader()
                csv_writer.writerow(params_dict)

# ---------------------- Data Analyser --------------------- #

import pandas as pd
import numpy as np

# animation
import matplotlib.pyplot as plt
import matplotlib.animation as animation
class DataAnalyser :
    def __init__(self, path_to_test_summary):
        self.df = pd.read_csv(path_to_test_summary+'.csv', sep = ',', header=0,index_colint=0) # take first line as headers


    def draw_speed_norm_distribution(self, test_id, save_animation = True):
        path_to_data = self.df.loc[[self.df['id_test'] == test_id,'path_to_data']]
        nb_parts = self.df.loc[[self.df['id_test'] == test_id,'total_number_of_particles']]
        self.df = pd.read_csv(path_to_data+'.csv', sep = ',', header=0,index_colint=0) # take first line as headers

        df['speed_norm'] = df.apply(compute_speed_norm, axis=1) # 1 for columns 

        # trying to create a group ... a group is {all particles at a time step}
        # In practice, we split the dataframe into several 
        row_count = len(DataFrame.index)
        lst = [df.iloc[i:i+nb_parts] for i in range(0,row_count-nb_parts+1,nb_parts)]
        number_of_frames = len(lst)
        
        def update_hist(num, lst):
            plt.cla()
            lst[num].hist['speed_norm']

        fig = plt.figure()
        fig.suptitle('{} : Speed distribution - {} iterations'.format(test_id,number_of_frames), fontsize=12)
        hist = lst[0].hist['speed_norm']
    
        animation = animation.FuncAnimation(fig, update_hist, number_of_frames, fargs=(lst, ) )
        if(save_animation):
            animation.save('{}_speed_norm_distribution_evolution.mp4'.format(test_id))
        plt.show()

    # --------------- utils ---------------- #
    def compute_speed_norm(row):
        return np.sqrt(row['vx']*row['vx']+row['vy']*row['vy']+row['vz']*row['vz'])





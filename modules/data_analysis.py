
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
        list_data=[]
        path = self.saving_directory + self.name_test
        if(not isfile(path + '.csv')):
            # does not already exist
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
import scipy.stats as ss

# animation
import matplotlib.pyplot as plt
import matplotlib.animation as animation
class DataAnalyser :
    def __init__(self, path_to_test_summary):
        self.df = pd.read_csv(path_to_test_summary+'.csv', sep = ',', header=0,index_col=0) # take first line as headers
        self.df.convert_dtypes().dtypes

    def draw_speed_norm_distribution(self, test_id, save_animation = True):
        path_to_data = self.df.loc[test_id].at['path_to_data'] # df.loc[5].at['B']
        nb_parts = self.df.loc[test_id].at['total_number_of_particles']
        df_data = pd.read_csv(path_to_data+'.csv', sep = ',', header=0, index_col=0) # take first line as headers

        # converting :
        df_data['vx'] = df_data['vx'].astype(float)
        df_data['vy'] = df_data['vy'].astype(float) 
        df_data['vz'] = df_data['vz'].astype(float) 

        df_data['speed_norm'] = df_data.apply(self.compute_speed_norm, axis=1) # 1 for columns 

        # trying to create a group ... a group is {all particles at a time step}
        # In practice, we split the dataframe into several 
        row_count = len(df_data.index)
        lst = [df_data.iloc[i:i+nb_parts] for i in range(0,row_count-nb_parts+1,nb_parts)]
        number_of_frames = len(lst)

        def update_hist(num, lst):
            #plt.cla() # to clear the current axis
            plt.clf() # to clear the current figure

            col = lst[num]['speed_norm']
            plt.hist(col, bins = 'auto', density=True)

            min_ = np.min(col)
            max_ = np.max(col)
            μ = np.mean(col) 
            σ = np.std(col)
            a = σ * np.sqrt(np.pi/(3.0*np.pi - 8.0)) # https://mathworld.wolfram.com/MaxwellDistribution.html
            m = 2.0*a*np.sqrt(2.0/np.pi)
            loc = μ - m

            X = np.linspace(min_, max_, 1000)
            Y = ss.maxwell.pdf(X, loc = loc, scale = a) 
            plt.plot(X,Y)

            fig.suptitle('{} : Speed distribution - iteration {}/{} - mean value {}'.format(test_id+1, num, number_of_frames,round(μ,1)), fontsize=12)

        fig = plt.figure(figsize=(10,10))

        
        col = lst[0]['speed_norm']
        plt.hist(col, bins = 'auto')

        min_ = np.min(col)
        max_ = np.max(col)
        μ = np.mean(col) # loc in the ss.maxwell function
        σ = np.std(col)
        a = σ * np.sqrt(np.pi/(3.0*np.pi - 8.0)) # https://mathworld.wolfram.com/MaxwellDistribution.html
        m = 2.0*a*np.sqrt(2.0/np.pi)
        loc = μ - m
        
        X = np.linspace(min_, max_, 1000)
        Y = ss.maxwell.pdf(X, loc = μ, scale = a) 
        plt.plot(X,Y)

        fig.suptitle('{} : Speed distribution - iteration {}/{} - mean value {}'.format(test_id, 1, number_of_frames,round(μ,1)), fontsize=12)

        interval = 200
        anim = animation.FuncAnimation(fig, update_hist, interval=interval, frames=number_of_frames, fargs=(lst, ), save_count=number_of_frames)
        if(save_animation):
            anim.save('{}_speed_norm_distribution_evolution.mp4'.format(test_id))
        else:
            plt.show()
    # --------------- utils ---------------- #
    def compute_speed_norm(self, row):
        return np.sqrt(row['vx']*row['vx']+row['vy']*row['vy']+row['vz']*row['vz'])





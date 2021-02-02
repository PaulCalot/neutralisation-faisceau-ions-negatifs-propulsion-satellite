
# what is a callback function : https://stackoverflow.com/questions/824234/what-is-a-callback-function

# NOTE : there is a warning appearing for reasons detailed here : 
# https://stackoverflow.com/questions/40659212/futurewarning-elementwise-comparison-failed-returning-scalar-but-in-the-futur
# there is nothing to do about it.

# csv usage : https://realpython.com/python-csv/#reading-csv-files-with-csv

import csv
from os.path import isfile, exists
import pandas as pd
import numpy as np
import scipy.stats as ss
from pathlib import Path, PurePath
from os import mkdir

# animation
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# my import
from main import get_maxwellian_params

class DataSaver :
    def __init__(self, particles, name_test, saving_directory):
        self.particles = particles
        self.name_test = name_test
        self.saving_directory = saving_directory

    def save_to_csv_(self, name, list_data, erase = True, use_saving_directory = True): # by default will erase what was already here
        
        path = self.saving_directory / name if use_saving_directory else Path(name)

        if(isfile(path) and not erase):
            mode = 'a'
        else :
            mode = 'w'
        with open(path, mode = mode) as csv_file:
                csv_writer = csv.writer(csv_file, delimiter=',')
                for data in list_data:
                    csv_writer.writerow(data)
    
    # def save_everything_to_multiple_csv(self, iter):
    #     # iter : iteration number
    #     list_data = [self.particles[0].get_headers()]
    #     for k in range(len(self.particles)):
    #         list_data.append(self.particles[k].to_list())
    #     self.save_to_csv_(self.name_test + "_iter{}".format(iter), list_data)

    def save_everything_to_one_csv(self, erase = False):
        list_data=[]
        path = self.saving_directory / (self.name_test)
        if(erase or not isfile(path)):
            # does not already exist
            list_data = [self.particles[0].get_headers()]
            
        for k in range(len(self.particles)):
            list_data.append(self.particles[k].to_list())
        
        self.save_to_csv_(path, list_data, erase = erase)

    # saving test params
    def save_test_params(self, tests_summary_file_name, params_dict, use_saving_directory = True, erase = False):
        path = self.saving_directory / (tests_summary_file_name) if use_saving_directory else Path(tests_summary_file_name)
        if(isfile(path) and not erase): # file already created
            df = pd.read_csv(path, sep = ',', header=0, index_col=0) # take first line as headers
            new_id = params_dict['id_test']
            count = 0
            L = [str(k) for k in list(df.index)]
            while(new_id in L):
                count+=1
                new_id = str(params_dict['id_test'])+'({})'.format(count)
            params_dict['id_test']=new_id
            df2 = pd.DataFrame(params_dict, index = [0], dtype = str)
            df2.set_index('id_test', inplace=True)
            df = df.append(df2)
        else:
            df = pd.DataFrame(params_dict, index = [0], dtype = str)
            #df.set_index('id_test', inplace=True)
            # df.astype({'id_test': str})
            # with open(path + '.csv', mode = 'w') as csv_file:
            #     fieldnames = params_dict.keys()
            #     csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            #     # first time, we add the keys
            #     csv_writer.writeheader()
            #     csv_writer.writerow(params_dict)
        df.to_csv(path) # should be updated

    def update_saved_params(self, tests_summary_file_name, params_dict, use_saving_directory = True):
        path = self.saving_directory / (tests_summary_file_name) if use_saving_directory else Path(tests_summary_file_name)
        if(exists(path)): # file already created
            df = pd.read_csv(path, sep = ',', header=0, index_col=0) # take first line as headers
            #df.set_index('id_test')
            df2 = pd.DataFrame(params_dict, index = [0], dtype = str)
            df.update(df2) # automatically update the row.
            df.to_csv(path) # should be updated
# ---------------------- Data Analyser ---------------------    
# To use only matplotlib: 
# https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/
# https://numpy.org/doc/stable/reference/generated/numpy.c_.html

# For pandas, several issues were met ... Name of the parameters, etc.


class DataAnalyser :
    def __init__(self, path_to_test_summary):   
        self.df = pd.read_csv(path_to_test_summary, sep = ',', header=0, index_col='id_test') # take first line as headers
        self.test_params = None
        self.current = None
        self.nb_parts = None
        self.number_of_frames = None
        self.test_id = None

    def load_test(self, test_id):
        self.test_id = test_id
        try : 
            self.test_params = self.df.loc[test_id]
        except KeyError:
            self.test_params = self.df.loc[int(test_id)]
        path_to_data = self.test_params.at['path_to_data'] # df.loc[5].at['B']
        self.path_to_saving = path_to_data.split('/')
        self.path_to_saving = '/'.join(self.path_to_saving[:len(self.path_to_saving)-1])+"/figures/"
        if(not exists(self.path_to_saving)): mkdir(self.path_to_saving)
        self.nb_parts = self.test_params.at['total_number_of_particles']
        df_data = pd.read_csv(path_to_data, sep = ',', header=0, index_col=0) # take first line as headers

        if(not 'speed_norm_squared' in df_data): # first time this test is loaded.
            print('First time test {} is loaded. Computing speed norm and saving it back.'.format(test_id))
            # converting data to int :
                # position
            df_data['x'] = df_data['x'].astype(float)
            df_data['y'] = df_data['y'].astype(float) 
            df_data['z'] = df_data['z'].astype(float)
                # speed
            df_data['vx'] = df_data['vx'].astype(float)
            df_data['vy'] = df_data['vy'].astype(float) 
            df_data['vz'] = df_data['vz'].astype(float) 

            # computing speed norm :
            df_data['speed_norm_squared'] = df_data.apply(self.compute_speed_norm_square, axis=1)
            df_data['speed_norm'] = df_data.apply(self.compute_speed_norm, axis=1) # 1 for columns 

            # updating the csv with the new columns
            df_data.to_csv(path_to_data) # should be updated

        # splitting in time step
        row_count = len(df_data.index)
        lst = [df_data.iloc[i:i+self.nb_parts] for i in range(0,row_count-self.nb_parts+1,self.nb_parts)]
        
        self.number_of_frames = len(lst)
        self.current = lst

        # --------------- Getter in csv files ----------------- #

    def get_param(self, key, test_id = None):
        # None can be used in case a test has already be defined
        id_row = test_id if test_id != None else self.test_id
        if(id_row == None):
            print('Please specify a test id.')
            return
        try : 
            return self.df.loc[id_row].at[key]
        except KeyError: 
            return self.df.loc[int(id_row)].at[key]

    # TODO
    def get_data(self, idx, time, key):
        pass
    # to get a line or a case for example, in the data.
    # it's doable even though can be a long getter

        # --------------- utils ---------------- #
    def compute_speed_norm(self, row):
        return np.sqrt(row['speed_norm_squared'])
                                                                      
    def compute_speed_norm_square(self, row):
        return (row['vx']*row['vx']+row['vy']*row['vy']+row['vz']*row['vz'])

        # -------------- drawing --------------- #
    def draw_hist_distribution(self, value_name = "", save_animation = True, plot_maxwellian = False, plot_gaussian = False, density = True, range = None, color = None):
        if(self.current == None):
            print("You have to load the test first using command : DataAnalyser.load_test(test_id)")
            return

        lst = self.current

        def update_hist(num, lst):
            plt.clf() # to clear the current figure # should not do that (that's what bugged after)

            col = lst[num][value_name]

            plt.hist(col, bins = 'auto', density=density, range = range, color = color)

            μ =np.mean(col)

            if(plot_maxwellian or plot_gaussian):
                min_ = np.min(col)
                max_ = np.max(col)
                X = np.linspace(min_, max_, 1000)
                if(plot_maxwellian):
                    loc, a = get_maxwellian_params(μ, np.std(col))
                    Y = ss.maxwell.pdf(X, loc = loc, scale = a)
                else :
                    Y = ss.norm.pdf(X, loc=μ, scale = np.std(col))
                plt.plot(X,Y, label = 'maxwellian pdf' if plot_maxwellian else 'gaussian pdf')

            ax.set_title('{} : {} distribution - iteration {}/{} - mean value {}'.format(self.test_id, value_name, num+1, self.number_of_frames,round(μ,1)), fontsize=12)
            plt.xlabel(value_name,fontsize=16)
            plt.ylabel("density",fontsize=16)

        fig = plt.figure(figsize=(10,10))
        
        col = lst[0][value_name]
        plt.hist(col, bins = 'auto', density=density, range = range, color = color)

        μ = np.mean(col)
        
        if(plot_maxwellian or plot_gaussian):
            min_ = np.min(col)
            max_ = np.max(col)
            X = np.linspace(min_, max_, 1000)
            if(plot_maxwellian):
                loc, a = get_maxwellian_params(μ, np.std(col))
                Y = ss.maxwell.pdf(X, loc = loc, scale = a)
            else :
                Y = ss.norm.pdf(X, loc=μ, scale = np.std(col))
            plt.plot(X,Y, label = 'maxwellian pdf' if plot_maxwellian else 'gaussian pdf')

        ax.set_title('{} : {} distribution - iteration {}/{} - mean value {}'.format(self.test_id, value_name, 1, self.number_of_frames,round(μ,1)), fontsize=12)
        plt.xlabel(value_name,fontsize=16)
        plt.ylabel("density",fontsize=16)
        interval = 200 # default value

        anim = animation.FuncAnimation(fig, update_hist, interval=interval, frames=self.number_of_frames, fargs=(lst, ), save_count=self.number_of_frames)
        if(save_animation):
            anim.save(self.path_to_saving+'{}_{}_hist_distribution.mp4'.format(self.test_id, value_name))
        else:
            plt.show()

    def draw_spatial_distribution(self, name = '' , save_animation = True, grid_size = 20, vmin = 0, vmax = 5e3):
        if(self.current == None):
            print("You have to load the test first using command : DataAnalyser.load_test(test_id)")
            return
    
        if(name == None):
            name = 'particles'

        lst = self.current
        def update_hist(num, lst, cb):
            plt.cla() # to clear the current figure
            df = lst[num]
            if(name != 'particles' and name != None):
                hexbin = ax.hexbin(df['x'], df['y'], df[name], gridsize = grid_size,  cmap='seismic', vmin = vmin, vmax = vmax)
            else :
                hexbin = ax.hexbin(df['x'], df['y'], None, gridsize = grid_size,  cmap='seismic', vmin = vmin, vmax = vmax)

            ax.set_title('{} : {} spatial distribution - iteration {}/{}'.format(self.test_id, name, num+1, self.number_of_frames), fontsize=12)
            
        fig, ax = plt.subplots(figsize=(10,10))
        df = lst[0]
        if(name != 'particles' and name != None):
            hexbin = ax.hexbin(df['x'], df['y'], df[name], gridsize = grid_size,  cmap='seismic', vmin = vmin, vmax = vmax)
        else :
            hexbin = ax.hexbin(df['x'], df['y'], None, gridsize = grid_size,  cmap='seismic', vmin = vmin, vmax = vmax)

        cb = fig.colorbar(hexbin, ax=ax)
        cb.set_label(name)

        ax.set_title(self.path_to_saving+'{} : {} spatial distribution - iteration {}/{}'.format(self.test_id,name, 1, self.number_of_frames), fontsize=12)
        

        interval = 200 # default value

        anim = animation.FuncAnimation(fig, update_hist, interval=interval, frames=self.number_of_frames, fargs=(lst, cb), save_count=self.number_of_frames)

        if(save_animation):
            anim.save(self.path_to_saving+'{}_{}_spatial_distribution_evolution.mp4'.format(self.test_id, name))
        else:
            plt.show()
      
    def draw_particles(self, save_animation=True):
        # useful : https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot
        # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.scatter.html 
        if(self.current == None):
            print("You have to load the test first using command : DataAnalyser.load_test(test_id)")
            return

        lst = self.current

        def update_hist(num, lst):
            #plt.cla() # to clear the current figure
            df = lst[num]
            # since we modifying scat we dont want to use plt.cla
            scat.set_offsets(np.c_[df['x'],df['y']])
            scat.set_array(df['speed_norm'])

            ax.set_title('{} : System evolution - iteration {}/{}'.format(self.test_id, num+1, self.number_of_frames), fontsize=15)

        fig, ax = plt.subplots(figsize=(10,10))
        df = lst[0]
        scat = ax.scatter(df['x'], df['y'], s=0.3, c = df['speed_norm'], cmap='seismic')

        ax.set_title('{} :  System evolution - iteration {}/{}'.format(self.test_id, 1, self.number_of_frames), fontsize=12)

        interval = 40 # default value

        anim = animation.FuncAnimation(fig, update_hist, interval=interval, frames=self.number_of_frames, fargs=(lst, ), save_count=self.number_of_frames)

        if(save_animation):
            anim.save(self.path_to_saving+'{}_system_evolution.mp4'.format(self.test_id))
        else:
            plt.show()

    # --------------------- Draw a frame ------------------------------ #

    def draw_hist_distribution_frame(self, which = 'last', value_name = "", save_frame = True, plot_maxwellian = False, plot_gaussian = False, density = True, range = None, color=None):
        if(self.current == None):
            print("You have to load the test first using command : DataAnalyser.load_test(test_id)")
            return

        lst = self.current

        frame = -1
        if(which == "first"):
            frame = 0
        elif(which == "last"):
            frame = len(lst)-1
        elif(type(which) == int):
            frame = which
        else:
            print("Please choose 'which' amongst : 'first', 'last' or any positive {}. Default is 'last'.".format(int))
            return

        df = lst[frame]
         
        fig, ax = plt.subplots(figsize=(10,10))
        
        col = df[value_name]
        plt.hist(col, bins = 'auto', density=density, range = range, color = color)

        μ = np.mean(col)
        std =  np.std(col)
        if(plot_maxwellian or plot_gaussian):
            min_ = np.min(col)
            max_ = np.max(col)
            X = np.linspace(min_, max_, 1000)
            if(plot_maxwellian):
                loc, a = get_maxwellian_params(μ, std)
                Y = ss.maxwell.pdf(X, loc = loc, scale = a)
            else :
                Y = ss.norm.pdf(X, loc=μ, scale = std)
            plt.plot(X,Y, label = 'maxwellian pdf' if plot_maxwellian else 'gaussian pdf')
        ax.set_title('{} : {} - it {}/{} - $\\mu$ : {} - $\\sigma$ : {}'.format(self.test_id, value_name, frame+1, self.number_of_frames,round(μ,1),round(std,2)), fontsize=15)
        plt.xlabel(value_name,fontsize=16)
        plt.ylabel("density",fontsize=16)
        plt.legend(loc='best')

        if(save_frame):
            plt.savefig(self.path_to_saving+'{}_{}_hist_distribution_it_{}.png'.format(self.test_id, value_name,frame+1), bbox_inches = 'tight', pad_inches = 0)
        else:
            plt.show()
        plt.close()

    def draw_spatial_distribution_frame(self, which = "last", name = '' , save_frame = True, grid_size = 20, vmin = 0, vmax = 5e3):
        # which : "first", "last" or any number between the 2.
        if(self.current == None):
            print("You have to load the test first using command : DataAnalyser.load_test(test_id)")
            return
    
        if(name == None):
            name = 'particles'

        lst = self.current
        

        fig, ax = plt.subplots(figsize=(10,10))
        frame = -1
        if(which == "first"):
            frame = 0
        elif(which == "last"):
            frame = len(lst)-1
        elif(type(which) == int):
            frame = which
        else:
            print("Please choose 'which' amongst : 'first', 'last' or any positive {}. Default is 'last'.".format(int))
            return
        df = lst[frame]

        if(name != 'particles' and name != None):
            hexbin = ax.hexbin(df['x'], df['y'], df[name], gridsize = grid_size,  cmap='seismic', vmin = vmin, vmax = vmax)
        else :
            hexbin = ax.hexbin(df['x'], df['y'], None, gridsize = grid_size,  cmap='seismic', vmin = vmin, vmax = vmax)

        cb = fig.colorbar(hexbin, ax=ax)
        cb.set_label(name)

        ax.set_title('{} : {} spatial distribution - iteration {}/{}'.format(self.test_id,name, frame +1, self.number_of_frames), fontsize=15)
        plt.legend(loc='best',fontsize=14)

        if(save_frame):
            plt.savefig(self.path_to_saving+'{}_{}_spatial_distribution_it_{}.png'.format(self.test_id, name, frame + 1), dpi = 300, bbox_inches = 'tight', pad_inches = 0)
        else:
            plt.show()
        plt.close()


    def draw_particles_frame(self, which = "last", save_frame=True, x_min = None, x_max = None, y_min = None, y_max = None):
        # useful : https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot
        # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.scatter.html 
        if(self.current == None):
            print("You have to load the test first using command : DataAnalyser.load_test(test_id)")
            return

        lst = self.current

        frame = -1
        if(which == "first"):
            frame = 0
        elif(which == "last"):
            frame = len(lst)-1
        elif(type(which) == int):
            frame = which
        else:
            print("Please choose 'which' amongst : 'first', 'last' or any positive {}. Default is 'last'.".format(int))
            return

        fig, ax = plt.subplots(figsize=(10,10))
        df = lst[frame]
        scat = ax.scatter(df['x'], df['y'], s=0.3, c = df['speed_norm'], cmap='seismic')
        ax.axis('equal')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        ax.set_title('{} :  System evolution - iteration {}/{}'.format(self.test_id, frame+1, self.number_of_frames), fontsize=15)

        plt.legend(loc='best',fontsize=14)

        if(save_frame):
            plt.savefig(self.path_to_saving+'{}_system_state_it_{}.png'.format(self.test_id, frame +1), dpi = 300, bbox_inches = 'tight', pad_inches = 0)
        else:
            plt.show()
        plt.close()

    # https://sciencing.com/maxwell-boltzmann-distribution-function-derivation-examples-13722756.html
    
    def savitzky_golay(self, y, window_size, order, deriv=0, rate=1):
        # https://fr.wikipedia.org/wiki/Algorithme_de_Savitzky-Golay
        
        import numpy as np
        from math import factorial

        try:
            window_size = np.abs(np.int(window_size))
            order = np.abs(np.int(order))
        except ValueError:
            raise ValueError("window_size and order have to be of type int")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        
        order_range = range(order+1)
        half_window = (window_size -1) // 2
        
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
        
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return np.convolve( m[::-1], y, mode='valid')

    def draw_Temperature_evolution(self, time_between_frames, tau_init, particles_mass, begin = 0, end = 1, save_frame = True):
        from scipy.optimize import least_squares

        if(self.current == None):
            print("You have to load the test first using command : DataAnalyser.load_test(test_id)")
            return

        lst = self.current
         
        fig = plt.figure(figsize=(10,10))
        Temp = []
        k = 1.380649e-23
        factor = particles_mass/(3*k) # m/(3k)
        begin_indx = int(begin*len(lst))
        end_indx = int(end*len(lst))
        number_of_frames_used = end_indx-begin_indx

        a = np.sqrt(np.mean(lst[0]['speed_norm_squared'])/3)
        eq_std = a*a*(3*np.pi-8)/(np.pi)
        print('Equilibrium variance = {} m2/s2'.format(eq_std))

        # time list 
        listTime = np.linspace(0,time_between_frames*number_of_frames_used,number_of_frames_used)
        
        # creating the list of temperatures
        for k in range(begin_indx,end_indx):
            df = lst[k]
            v = df['speed_norm']
            v_mean = np.mean(v)
            w = v - v_mean
            variance_t = np.mean(w*w)
            # 
            Temp.append(variance_t)
        Temp = np.array(Temp)
        T0 = Temp[0]
        print('T_0 = {} K'.format(T0))
        # minimization problem

        def f(X_, Temp,listTime, T0): # X = [tau]
            tau = X_[0]
            f_ = lambda T, t: np.abs(T - (T0-eq_std)*np.exp(-t/tau)-eq_std)
            total = 0
            for k in range(len(Temp)):
                T = Temp[k]
                t = listTime[k]
                f_value = f_(T,t)
                total+=f_value*f_value
            return total
        #Temp_smooth = self.savitzky_golay(np.array(Temp), window_size = 13, order=3)
        results = least_squares(f, np.array([tau_init]), bounds = ([1e-3*tau_init],[1e3*tau_init]), args = (Temp,listTime,T0)).x
        tau = results[0]

        def get_Temp(Time):
            return (T0-eq_std)*np.exp(-Time/tau)+eq_std

        fig, ax = plt.subplots(figsize=(15,10))
        ax.set_title("évolution de la variance de la norme de la vitesse - $\sigma_e$ = {} m/s ; $\\tau$ = {} s".format("{:e}".format(np.sqrt(eq_std)),"{:e}".format(tau)))
        plt.xlabel("temps (s)",fontsize=16)
        plt.ylabel("variance de la vitesse ($m^2/s^2$)",fontsize=16)
        plt.plot(listTime, get_Temp(listTime), label = '$\sigma^2(t) = (\sigma^2(0)-\sigma^2_e)exp(-t/\\tau)+\sigma^2_e$')
        plt.plot(listTime,Temp, label = 'Simulation')
        #plt.plot(listTime,Temp_smooth)
        plt.legend(loc='best',fontsize=14)
        if(save_frame):
            plt.savefig(self.path_to_saving+'{}_temperature_evolution.png'.format(self.test_id), bbox_inches = 'tight', pad_inches = 0)
        else:
            plt.show()
        plt.close()

    def draw_Temperature_evolution_per_axis(self, time_between_frames, tau_init, particles_mass, begin = 0, end = 1, save_frame = True):
        from scipy.optimize import least_squares

        if(self.current == None):
            print("You have to load the test first using command : DataAnalyser.load_test(test_id)")
            return

        lst = self.current
         
        k = 1.380649e-23
        factor = particles_mass/(3*k) # m/(3k)
        begin_indx = int(begin*len(lst))
        end_indx = int(end*len(lst))
        number_of_frames_used = end_indx-begin_indx

        # teq 
        a = np.sqrt(np.mean(lst[0]['speed_norm_squared'])/3)
        T_eq = a*a*particles_mass/k
        print('Teq = {} K'.format(T_eq))

        # time list 
        listTime = np.linspace(0,time_between_frames*number_of_frames_used,number_of_frames_used)
        TX, TY, TZ = [], [], []
        # creating the list of temperatures
        for k in range(begin_indx,end_indx):
            df = lst[k]
            vx = np.std(df['vx'])
            vy = np.std(df['vy'])
            vz = np.std(df['vz'])
            TX.append(factor*vx*vx)
            TY.append(factor*vy*vy)
            TZ.append(factor*vz*vz)
        T0X,T0Y,T0Z = TX[0],TY[0],TZ[0]
        T0 = [T0X,T0Y,T0Z]
        Temp = [TX,TY,TZ]
        # minimization problem
        def f(X_, Temp,listTime, T0): # X = [tau]
            tau = X_[0]
            f_ = lambda T, t: np.abs(T - (T0-T_eq)*np.exp(-t/tau)-T_eq)
            total = 0
            for k in range(len(Temp)):
                T = Temp[k]
                t = listTime[k]
                f_value = f_(T,t)
                total+=f_value*f_value
            return total
        #Temp_smooth = self.savitzky_golay(np.array(Temp), window_size = 13, order=3)
        results_X = least_squares(f, np.array([tau_init]), bounds = ([1e-3*tau_init],[1e3*tau_init]), args = (Temp,listTime,TX)).x        
        results_Y = least_squares(f, np.array([tau_init]), bounds = ([1e-3*tau_init],[1e3*tau_init]), args = (Temp,listTime,TY)).x
        results_Z = least_squares(f, np.array([tau_init]), bounds = ([1e-3*tau_init],[1e3*tau_init]), args = (Temp,listTime,TZ)).x
        results = [results_X, results_Y, results_Z]
        tau_X = results_X[0]
        tau_Y = results_Y[0]
        tau_Z = results_Z[0]
        tau = [tau_X,tau_Y, tau_Z]

        def get_Temp(Time, T0, tau):
            return (T0-T_eq)*np.exp(-Time/tau)+T_eq
        axis = ['x','y','z']
        for k in range(3):
            fig, ax = plt.subplots(figsize=(15,10))
            ax.set_title("Thermal evolution along axis {} - $\\tau$ = {} s".format(axis[k], "{:e}".format(tau[k])))
            plt.xlabel("time (s)",fontsize=16)
            plt.ylabel("Temperature (K)",fontsize=16)
            plt.plot(listTime, get_Temp(listTime ,T0[k], tau[k]), label = '$T^2(t) = (T^2(0)-T^2_e)exp(-t/\\tau)+T^2_e$')
            plt.plot(listTime, Temp[k], label = 'Simulation')
            #plt.plot(listTime,Temp_smooth)
            plt.legend(loc='best',fontsize=14)
            if(save_frame):
                plt.savefig(self.path_to_saving+'{}_temperature_evolution_{}.png'.format(self.test_id, axis[k]), bbox_inches = 'tight', pad_inches = 0)
            else:
                plt.show()
            plt.close()

def merge_tests_summary(names, output_name):
    mode = "r"
    L = []
    for k, name in enumerate(names) : 
        with open(name + '.csv', mode = mode) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')            
            for i, row in enumerate(csv_reader):
                if(i>0 or k==0): # we do not take the first one which is the raw with the header
                    L.append(row)
    
    with open(output_name + '.csv', mode = 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')            
        for row in L:
            csv_writer.writerow(row)
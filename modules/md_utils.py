import numpy as np
import os
from numpy.core.fromnumeric import var
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def complementaire(i):
    if i+1<10:
        return "000"+str(i+1)
    elif i+1<100:
        return "00"+str(i+1)
    elif i+1<1000:
        return "0"+str(i+1)
    else :
        return str(i+1)
    
def lenght_file(path):
    with open(path, "r") as file:
        ind=0
        line = file.readline()
        while line!="":
            line = file.readline()
            ind+=1
        file.close()
        return ind

def lenght_headers(path):
    with open(path, "r") as file:
        file = open(path, "r")
        ind=0
        line = file.readline()
        while line[0]=="#" or line[0]=="%":
            line = file.readline()
            ind+=1
        file.close()
        return ind


def scan_directory(dir_path, file_extension):
        number_of_files = 0
        names = []
        for file in os.listdir(dir_path):
            if file.endswith(file_extension):
                number_of_files+=1
                names.append(file)
        return names, number_of_files


def get_impact_time(df):
    results = np.where(df['imp?']==1)[0] # 1st time imp? = 1
    if len(results)>0: return results[0]
    else: return None

# ------------------------------------ Plotting ----------------------------------- #

def plot_crystal(ax, dataframe, radii, colors, dir1, dir2, alpha = 0.1):
    for index, row in dataframe.iterrows():
        ax.add_patch(plt.Circle((row[dir1], row[dir2]), radii[row['type']], \
                                color = colors[row['type']], alpha = alpha))

def plot_2d_simu(type_ion, df_ion, df_crystal, radii, colors, start = True):
    
    xmin, xmax = min(np.min(df_crystal['x']),np.min(df_ion['x'])), max(np.max(df_crystal['x']),np.max(df_ion['x']))
    ymin, ymax = min(np.min(df_crystal['y']),np.min(df_ion['y'])), max(np.max(df_crystal['y']),np.max(df_ion['y']))
    zmin, zmax = min(np.min(df_crystal['z']),np.min(df_ion['z'])), max(np.max(df_crystal['z']),np.max(df_ion['z']))

    fig, axes = plt.subplots(1,2, figsize = (20,10))
    axes[0].axis('equal')
    axes[0].set(xlim=(xmin, xmax), ylim = (zmin, zmax))
    plot_crystal(axes[0], df_crystal, radii, colors, 'x', 'z', alpha = 0.1)
    axes[0].set_xlabel(r'x ($\AA$)', fontsize = 14)
    axes[0].set_ylabel(r'z ($\AA$)', fontsize = 14)

    axes[1].axis('equal')
    axes[1].set(xlim=(ymin, ymax), ylim = (zmin, zmax))
    plot_crystal(axes[1], df_crystal, radii, colors, 'y', 'z', alpha = 0.1)
    axes[1].set_xlabel(r'y ($\AA$)', fontsize = 14)
    axes[1].set_ylabel(r'z ($\AA$)', fontsize = 14)
    
    if(start):
        axes[0].add_patch(plt.Circle((df_ion['x'][0],df_ion['z'][0]), radii[type_ion], color = colors[type_ion]))
        axes[0].plot(df_ion['x'], df_ion['z'], color = colors[type_ion], linewidth = 3)

        axes[1].add_patch(plt.Circle((df_ion['y'][0],df_ion['z'][0]), radii[type_ion], color = colors[type_ion]))
        axes[1].plot(df_ion['y'], df_ion['z'], color = colors[type_ion], linewidth = 3)
    else:
        l = len(df_ion['x'])
        axes[0].add_artist(plt.Circle((df_ion['x'][l-1],df_ion['z'][l-1]), radii[type_ion], color = colors[type_ion]))
        axes[1].add_artist(plt.Circle((df_ion['y'][l-1],df_ion['z'][l-1]), radii[type_ion], color = colors[type_ion]))
    
    impact_time = get_impact_time(df_ion)
    if(impact_time is not None):
        # there has been an impact
        row_impact = df_ion.iloc[impact_time]
        print('Impact time : {} ps.'.format(row_impact['Integration time']))
        axes[0].plot(row_impact['x'],row_impact['z'], 'x', color = 'r', markersize = 20)
        axes[1].plot(row_impact['y'],row_impact['z'], 'x', color = 'r', markersize = 20)

    plt.show()


# ------------------------------------------- Reading results folder ------------------------------- #

class process_results:
    def __init__(self, path_to_folder):
        self.path = path_to_folder # path to folder
        self.ion = process_ion_folder(path_to_folder/'ion')
        self.crystal = process_crystal_folder(path_to_folder/'cfg')
        self.log = process_log_file(path_to_folder/'log')

# ------------------------------------------- Reading ion folder ----------------------------------- #

class process_ion_folder :
    # will process everything in the ion folder (every file with .ion)
    file_fields = ['st', 't', 'ke/E_i', 'ipe', 'epe', 'd.x', 'd.y', 'd.z', '|d|', '#b', 'imp?']
    final_fields = ['st', 't', 'ke/E_i', 'ipe', 'epe', 'x', 'y', 'z', '|d|', 'imp?']
    
    def __init__(self, path):
        self.path = path
        self.dataframes = self.get_dict_df()        
        self.list_angles = self.extract_angles()

    def get_dict_df(self):
        self.names, self.simulation_number = scan_directory(self.path, '.ion')
        df_dict = {}
        for name in self.names :
            df_dict[name] = self.get_dataframe_(self.path/name)
        return df_dict

    def get_dataframe_(self, path):
        # wait so the ion file can be written before executing this script.
        l = lenght_headers(path)
        df=pd.read_csv(path, header=None, skiprows=l, names=self.file_fields, sep=" ")

        xini,yini,zini = self.coord_initiales_ion_(path)

        df['x']=xini+df['d.x']
        df['y']=yini+df['d.y']
        df['z']=zini+df['d.z']

        df.drop(['#b','d.x', 'd.y','d.z'], axis='columns', inplace=True)
        df.set_index('st') # step time as index
        
        return df

    def extract_angles(self):
        list_angles = []
        for k, df in self.dataframes.items():
            dx= df['x'].values[-1]- df['x'].values[-2]
            dy= df['y'].values[-1]- df['y'].values[-2]
            dz= df['z'].values[-1]- df['z'].values[-2]
            dr=np.sqrt(dx**2+dy**2+dz**2)

            phi_sortie = np.nan
            theta_sortie=np.arccos(dz/dr)*180/np.pi
            if(dx != 0):
                phi_sortie=np.arctan(dy/dx)*180/np.pi
            list_angles.append([theta_sortie, phi_sortie])
        return list_angles
    
    def coord_initiales_ion_(self, path):       
        with open(path, "r") as file:
            found = False
            while not found:
                line = file.readline()

                if(line == ""): # https://docs.python.org/3.6/tutorial/inputoutput.html#methods-of-file-objects
                    print('End of file reached. Cound not find initial ion coordinates.')
                    return

                split = line.split(':')
                if split[0]=="# Initial COM position":
                    file.close()
                    found = True
                    xini,yini,zini = split[1].strip().split(' ')
                    return float(xini), float(yini), float(zini)
                
# ------------------------------------------  Reading crystal folder (cfg) ------------------------------- #

class process_crystal_folder :
    file_fields = ['id','type','x',\
        'y','z','vx','vy','vz','fix ?']

    def __init__(self, path):
        self.path = path
        self.dataframes = self.get_dict_df()

    def get_dict_df(self):
        self.names, self.simulation_number = scan_directory(self.path, '.cfg')
        df_dict = {}
        for name in self.names :
            df_dict[name] = self.get_df_crystal_(self.path/name)
        return df_dict

    def get_df_crystal_(self, path):
        df_crystal=pd.read_csv(path, header=None, skiprows=lenght_headers(path), names=self.file_fields, sep="\t")
        return df_crystal

    def plot_crystal(self, ax, radii, colors, dir1, dir2, idx_or_name = 0):
        if(type(idx_or_name) == int):
            plot_crystal(ax, self.dataframes[self.names[idx_or_name]], radii, colors, dir1, dir2, alpha = 0.1)
        elif(idx_or_name in self.names):
            plot_crystal(ax, self.dataframes[idx_or_name], radii, colors, dir1, dir2, alpha = 0.1)
        else:
            print("Name / index not recognized. Should be in {} / lenght [0, {}].".format(self.names, len(self.names)))
# ------------------------------------------- Reading log file ------------------------------- #


class process_log_file:
    
    def __init__(self, path):
        self.path = path
        self.integration_data, self.other_data, self.nb_traj  = self.get_data()    
    
    def get_data(self):
        int_data = []
        other_data = []
        with open(self.path, "r") as file:
            line = file.readline()
            if(line[0] == "#"):
                other_data.append(line)
            
            idx_traj = -1
            traj_data = []
            while(line!=""):
                line = file.readline()
                if(len(line)==0):
                    continue

                if(line[0] == "#"):
                    if(line[2:12] == 'Trajectory'):
                        if(idx_traj>-1):
                            int_data.append(traj_data)
                            traj_data = []
                        idx_traj+=1
                    other_data.append(line)
                else:
                    values_str = line.strip().split('\t')
                    traj_data.append([float(value) for value in values_str])
        headers = ['t [#]', 'time [ps]', 'KE [eV]', 'PE [eV]', 'E [eV]', 'Edrift [eV]', '%Edrift [%]', 'T [K]', 'nClst [#]', 'nAicl? [0/1]', 'Plat [GPa]']
        df = [pd.DataFrame(int_data[k], columns = headers) for k in range(idx_traj)]
        return df, other_data, idx_traj

    def plot_instantaneous_temperature_K(self, idx_traj = 0):
        self.plot_evolution(idx_traj, 'T [K]')
        # plt.plot(self.integration_data[idx_traj]['time [ps]'], self.integration_data[idx_traj]['T [K]'])

    def plot_evolution(self, idx_traj, y):
        # y = f(t)
        self.plot(idx_traj, 'time [ps]', y)
        
    def plot(self, idx_traj, x, y):
        # plot y as a function of x
        # where x and y should be either indexes or column labels.
        plt.plot(self.integration_data[idx_traj][x], self.integration_data[idx_traj][y])
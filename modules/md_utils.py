import numpy as np
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def coord_initiales_ion(path):
    lenght_headers = longueur_file(path)
    file = open(path, "r")
    
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

    # indice_2_points=0
    # indice_1er_espace=0
    # indice_2eme_espace=0
    # for i in range(len(line)):
    #     if indice_2_points==0 and line[i]==':':
    #         indice_2_points=i
    #     elif indice_1er_espace==0 and indice_2_points!=0 and i>indice_2_points+1 and line[i]==' ':
    #         indice_1er_espace=i
    #     elif indice_2eme_espace==0 and indice_1er_espace!=0 and line[i]==' ':
    #         print(i)
    #         indice_2eme_espace=i
            
    # xini=float(line[indice_2_points+2: indice_1er_espace])
    # yini=float(line[indice_1er_espace+1: indice_2eme_espace])
    # zini=float(line[indice_2eme_espace+1:])
    
    # return xini,yini,zini

def complementaire(i):
    if i+1<10:
        return "000"+str(i+1)
    elif i+1<100:
        return "00"+str(i+1)
    elif i+1<1000:
        return "0"+str(i+1)
    else :
        return str(i+1)
    
def longueur_file(chemin):
    file = open(chemin, "r")
    ind=0
    line = file.readline()
    while line!="":
        line = file.readline()
        ind+=1
    file.close()
    return ind

def longueur_intro(chemin):
    file = open(chemin, "r")
    ind=0
    line = file.readline()
    while line[0]=="#" or line[0]=="%":
        line = file.readline()
        ind+=1
    file.close()
    return ind


# ---------------------------------------- Extract data ---------------------------------- #
def extract_data_df(nb_steps, path):
    l = longueur_intro(path)
    df=pd.read_csv(path, header=None, skiprows=l, names=['Time step','Integration time','KE/TE_i','internal PE','external PE','d.x','d.y','d.z','norm(d)','#b','imp?'], sep=" ")
    while(nb_steps+1 < len(df.values)): # because the file is not necessarily written yet even though the previous cell has been completed.
        df=pd.read_csv(path, header=None, skiprows=l, names=['Time step','Integration time','KE/TE_i','internal PE','external PE','d.x','d.y','d.z','norm(d)','#b','imp?'], sep=" ")
   
    xini,yini,zini = coord_initiales_ion(path)
    df['x']=xini+df['d.x']
    df['y']=yini+df['d.y']
    df['z']=zini+df['d.z']

    df.drop(['#b','d.x', 'd.y','d.z'], axis='columns', inplace=True)
    df.set_index('Time step')
    return df

def extract_angles(df):
    dx=df['x'].values[-1]-df['x'].values[-2]
    dy=df['y'].values[-1]-df['y'].values[-2]
    dz=df['z'].values[-1]-df['z'].values[-2]
    dr=np.sqrt(dx**2+dy**2+dz**2)

    phi_sortie = np.nan
    theta_sortie=np.arccos(dz/dr)*180/np.pi
    if(dx != 0):
        phi_sortie=np.arctan(dy/dx)*180/np.pi
    else :
        print('Ion did not leave.')
    return theta_sortie, phi_sortie

def get_df_crystal(path):
    cols = ['id','type','x',\
       'y','z','vx','vy','vz','fix ?']
    df_crystal=pd.read_csv(path, header=None, skiprows=longueur_intro(path), names=cols, sep="\t")
    return df_crystal

def get_impact_time(df):
    results = np.where(df['imp?']==1)[0] # 1st time imp? = 1
    if len(results)>0: return results[0]
    else: return None
# ------------------------------------ Plotting ----------------------------------- #
def plot_crystal(ax, dataframe, radii, colors, dir1, dir2, alpha = 0.1):
    for index, row in dataframe.iterrows():
        ax.add_patch(plt.Circle((row[dir1], row[dir2]), radii[row['type']], \
                                color = colors[row['type']], alpha = alpha))

def plot_2D_simu(type_ion, df_ion, df_crystal, radii, colors, start = True):
    
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

# TODO - not sure it is useful though
# ------------------------------------------- Reading results folder ------------------------------- #

# ------------------------------------------- Reading ion folder ----------------------------------- #

# ------------------------------------------  Reading crystal folder (cfg) ------------------------------- #

# ------------------------------------------- Reading log file ------------------------------- #

class process_log_file:
    
    def __init__(self, path):
        self.path = path
        self.integration_data, self.other_data  = self.get_data(path)    
    
    def get_data(self, path):
        int_data = []
        other_data = []
        with open(path, "r") as file:
            line = file.readline()
            if(line[0] == "#"):
                other_data.append(line)
            
            idx_traj = -1
            traj_data = []
            while(line!=""):
                line = file.readline()
                
                if(len(line)>0):
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
        headers = ['t [#]', 'time [ps]', 'KE [eV]', 'PE [eV]', 'E [Ev]', 'Edrift [eV]', '%Edrift [%]', 'T [K]', 'nClst [#]', 'nAicl? [0/1]', 'Plat [GPa]']
        df = [pd.DataFrame(int_data[k], columns = headers) for k in range(idx_traj)]
        return df, other_data

    def plot_instantaneous_temperature_K(self, idx_traj = 0):
        plt.plot(self.integration_data[idx_traj]['time [ps]'], self.integration_data[idx_traj]['T [K]'])
        plt.show()
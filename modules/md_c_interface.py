import os
from pathlib import Path

def clean_and_remake(src_path, paths_to_h):
    #os.chdir(src_path)
    os.system("make -C {}/ clean".format(src_path))
    cmd = ''
    if(paths_to_h is not None):
        for path in paths_to_h:
            cmd += "-I{}".format(path)
    os.system("make {} -C {}/ md2".format(cmd, src_path))
    
def init_simulation(src_path, name):
    test_path = src_path/name
    # creates a folder in the src folder, named *name, for the new simulation
    if os.path.exists(test_path):
        os.system('rm -r {}'.format(test_path))
    os.makedirs(test_path)
    os.makedirs(test_path/'ion')
    os.makedirs(test_path/'cfg')
    os.makedirs(test_path/'clu')
    
    cmd1 ="cp -R {}/ {}".format(src_path/'cfg', test_path) # copying directory (-R) cp -R source target
    cmd2 = "cp {} {}".format(src_path/'md2', test_path/'md2')
    os.system(cmd1)
    os.system(cmd2)
    return test_path

def make_command(params, flags):
    possible_keys = ['-ion','-ionE','-Tset','-ionT','-ionP','-tau','-n','-dt','-i1']
    possible_flags = ['+dtv']
    cmd = ''
    
    for (key, value) in params.items():
        if(key not in possible_keys):
            print('Parameter should be in {}. But {} is not.'.format(possible_keys, key))
            break
        cmd += '{} {} '.format(key,      value)
        
    for flag in flags:
        if(flag not in possible_flags):
            print('Flag should be in {} but {} is not'.format(possible_flags, flag))
            break
        cmd += '{} '.format(flag)
    command = "./md2 -oc cfg/####.cfg {} > log &".format(cmd)
    return command

def launch_simulation(name, params, flags, paths_to_h = None):
    # name : name of the simulation (under the source folder of the md code)
    # params : a dictionnary containing as keys the items in the command to launch a simu in the C code
    
    current_dir = Path.cwd()
    src_path = current_dir/'md2_sources' # path to src folder

    # Preparing config - hypothesis : this notebook is in the same directory that the source MD code.

    clean_and_remake(src_path, paths_to_h)

    # all is saved in a new file
    test_path = init_simulation(src_path, name)
    
    # issuing command
    command = make_command(params, flags)
    print(command)
    os.chdir(test_path)
    os.system(command)
    os.chdir(current_dir)
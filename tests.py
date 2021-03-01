# local import
from main import simulate, processing_only

# python import
import sys, getopt
from pathlib import Path
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"ht:a:s:p:",['type=', 'archi=', 'sim=', 'proc=']) # h for help
    except getopt.GetoptError:
        print('Usage : ')
        print('main.py -t <type> -a <path_to_cfg_architecture> -s <path_to_cfg_simulation> -p <path_to_cfg_processing>')   
        sys.exit(2)

    if(len(opts)==0):
        print('Usage : ')
        print('main.py -t <type> -a <path_to_cfg_architecture> -s <path_to_cfg_simulation> -p <path_to_cfg_processing>')   
        sys.exit()

    if(opts[0][0] == '-h'):
        print('Usage : ')
        print('main.py -t <type> -a <path_to_cfg_architecture> -s <path_to_cfg_simulation> -p <path_to_cfg_processing>')   
        sys.exit()

    elif opts[0][0] in ('-t','--type'):
        # resolve allows to go in absolute path (and note relative anymore)
        mod_path = Path(__file__).resolve().parents[0]
        #print(opts)
        if(opts[1][0] in ('-a','--archi') and \
            opts[2][0] in ('-s','--sim') and \
                opts[3][0] in ('-p','--proc')):
            system = (mod_path/Path(opts[1][1])).resolve()
            simu = (mod_path/Path(opts[2][1])).resolve()
            proc = (mod_path/Path(opts[3][1])).resolve()
            #print(system, simu, proc)
            if(opts[0][1] == 'simulation'):
                simulate(system, simu, proc)
            elif(opts[0][1]  == 'processing'):
                processing_only(system, simu, proc)
        else :
            print('Please precise path to cfg files.')
            print('Usage : ')
            print('main.py -t <type> -a <path_to_cfg_architecture> -s <path_to_cfg_simulation> -p <path_to_cfg_processing>')   
            sys.exit()
    else :
        pass

if __name__ == "__main__":
    main(sys.argv[1:]) # because the first arg is always the name of the file
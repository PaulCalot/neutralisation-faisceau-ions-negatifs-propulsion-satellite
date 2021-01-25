# ------------------------------ Reading cfg files ---------------- #

def str_to_bool(s):
    """Convert Strings like "True" and "False" to booleans.
    Args:
        s (String): A string that should be "True" or "False"
    Raises:
        ValueError: If String has a wrong value
    Returns:
        (bool): Converted bool.
    """
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError
            
def read_args_multiple_types(input):
    # we say possible args are : int, bool, float, none and str
    if(input=='False' or input == 'True'):
        return str_to_bool(input)
    elif(input=='None'):
        return None
    try :
        integer = int(input)
        if(float(input)%integer==0):
            # it really is a integer
            return integer
        else:
            return float(intput)
    except : 
        return str(input)

def get_next_valid_idx(list, k):
    j = k+1
    while( j < len(list) and list[j]<list[k]):
        j+=1
    return j

class Options(object):

    def __init__(self, cfg):
        '''Defines the cfg file'''
        self.cfg = cfg

# # TODO: make use_dict a list so we can choose weather to send a list or a dict each time.
# def get_list_args(args, idx_args, use_dict = False, name_args = None, verbose = True):
#     """ This function takes a list of arguments and returns a list of lists or dicts where arguments are splitted as given in idx_args. 
#     Args:
#         args (list): list of arguments of any given type that are to be grouped together according to the idx_arg list.
#         idx_args (list of int): list of int which indicated how many args are for the current group.
#         use_dict (bool, optional): weather to use dict or not. Defaults to False.
#         name_args (list, optional): The name of the args if we use a dict (if we use a list, then the args are passed in the given order and no name is necessary). Defaults to None.
    
#     Returns:
#         (list) : a list of dict or list where each one is grouped up according to idx_args. 
#     """
#     n = len(idx_args) # the number of groups of args to make.
    
#     if(verbose):
#         print('Creating list of arguments ...')
#         print('')
#         print("Arguments : {}".format(args))
#         print("Group of args indexes : {}".format(idx_args))
#         print("Use dictionnary : {}".format(use_dict))
#         if(name_args != None):
#             print("Name of the arguments : {}".format(name_args))
#         print("Considering the data, {} groups of args will be made.".format(n))
    
#     layers_args = []
#     count = 0
#     count_dict = 0
#     for k in range(n):
#         if(verbose) : print('')
#         if(idx_args[k] == 0): # no args or autoset (it depends on you)
#             if(verbose) : print("Group #{} has no arg associated.".format(k))
#             arg = None
#         else:
#             if(verbose):
#                 print("Group #{} has {} args associated.".format(k,idx_args[k]) if idx_args[k] > 1 else "group #{} has {} arg associated.".format(k,idx_args[k]))

#             if(use_dict[k] if type(use_dict) != bool else use_dict):
#                 arg = dict()
#                 for j in range(count,count+idx_args[k]):
#                     arg[name_args[count_dict+j-count]] = args[j]
#                     if(verbose) : 
#                         print("New key : {}".format(name_args[count_dict+j-count]))
#                         print('Current args state : {}'.format(arg))
#                 count_dict += idx_args[k]
#             else :
#                 arg = []
#                 for j in range(count,count+idx_args[k]):
#                     arg.append(args[j])
#                     if(verbose) : print('Current args state : {}'.format(arg))
#             count += idx_args[k]
#             if(verbose) : print("Current arguments count : {}".format(count))
#         layers_args.append(arg)

#     if(verbose):
#         print('')
#         print("Output : {}".format(layers_args))
        
#     return layers_args
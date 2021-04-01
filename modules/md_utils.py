import numpy as np
import os


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
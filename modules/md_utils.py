import numpy as np
import os


def coord_initiales_ion(chemin):
    file = open(chemin, "r")
    
    test=True
    while test: #ou 8 ?
        line = file.readline()
        extrait=line[10:13]
        if extrait=="COM":
            test=False
    file.close()
    indice_2_points=0
    indice_1er_espace=0
    indice_2eme_espace=0
    for i in range(len(line)):
        if indice_2_points==0 and line[i]==':':
            indice_2_points=i
        elif indice_1er_espace==0 and indice_2_points!=0 and i>indice_2_points+1 and line[i]==' ':
            indice_1er_espace=i
        elif indice_2eme_espace==0 and indice_1er_espace!=0 and line[i]==' ':
            print(i)
            indice_2eme_espace=i
            
    xini=float(line[indice_2_points+2: indice_1er_espace])
    yini=float(line[indice_1er_espace+1: indice_2eme_espace])
    zini=float(line[indice_2eme_espace+1:])
    
    return xini,yini,zini

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
    while line[0]!="0":
        line = file.readline()
        ind+=1
    file.close()
    return ind-1
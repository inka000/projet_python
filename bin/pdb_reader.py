'''pdb_reader module
pdb_reader module allows to read a pdb file and to get the ATOM lines in a list of Atome class instances
'''
#!/usr/bin/python3

import re
import sys
import os
from classes import *




def getAtome(line) :
    '''
    Gets atoms informations from a string and creates an Atome instance
    '''
    return (Atome(line[13:16].split(" ")[0], str(line[21]), float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20], int(line[22:26]), int(line[6:11])))


def readPDB(filename):
    '''
    Reads the provided file in the data directory, uses only the first model
    If the file provided does not exist, it informs the user
    '''
    f=open("../data/"+filename,"r")
    atomes=[]
    re_atomes=re.compile("^ATOM")
    re_fin_modele=re.compile("^ENDMDL")
    for line in f:
        if re_atomes.search(line):
            atomes.append(getAtome(line))
        elif re_fin_modele.search(line):
            break
    f.close()
    return atomes



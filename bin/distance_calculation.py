'''distance_calculation module
distance_calculation module allows to calculate distance between two Atomes instances
'''

#!/usr/bin/python3

from math import sqrt
from pdb_reader import *
from classes import *



def distance(atom1, atom2):
    '''
    Calculates and returns the distance between two Atome instances based on coordinates, using sqrt function from math module
    ''' 
    d=sqrt((atom1.xpos-atom2.xpos)**2+(atom1.ypos-atom2.ypos)**2+(atom1.zpos-atom2.zpos)**2)
    return d

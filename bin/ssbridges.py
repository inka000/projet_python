'''ssbridges module
ssbridges module calculates disulfide bridges in a protein calling ssbridges_call function
'''


#!/usr/bin/python3

from recup import *
from verif_distance_diS import *

def ssbridge_call(dico, atomes, dico_dist, intra_inter_choice):
    '''
    ssbridge_call function calculates disulfide bridges in a protein. It takes in input a list of Atome instances which constitute the protein. This list have to be created by the pdb_reader module.
    Usefull information are indicated in  :
        - dico is a dictionary containing atoms of interest that may participate to a specified interaction type
        - dico_dist is a dictionary containing distance cut-off for each interaction type
    The module's function first gets atoms of interest with recup_S function of recup module.
    Then, it tries to calculate distances and chooses distances under 5.3 Angstrom with distance_file function from verif_distance module. 
    If there is not any atom of interest, an exception is raised and the interactions list is empty
    If there is not any distance according to cut-off distances, the list is also empty
    '''
    liste_S=recup_S(dico, atomes) 
    try :
        liste_of_DiS=distance_file(liste_S, dico_dist, "SSBOND", intra_inter_choice) 
        liste_of_interactions=(liste_of_DiS)

    except :
        liste_of_interactions=[]
        print("There is not any SSBOND")

    return liste_of_interactions
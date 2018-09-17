'''ionic module
ionic module calculates ionic interactions in a protein calling ionic_call function
'''

#!/usr/bin/python3

from recup import *
from verif_distance_diS import *

def ionic_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice):
    '''
    ionic_call function calculates ionic interactions in a protein. It takes in input a list of Atome instances which constitute the protein. This list have to be created by the pdb_reader module.
    Usefull information are indicated :
        - dico is a dictionary containing residues of interest that may participate to a specified interaction type
        - dico_atoms is a dictionary containing atoms of interest that may participate to a specified interaction type
        - dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    The module's function first gets atoms of interest with recup_ionics function of recup module.
    Then, it tries to calculate distances and chooses distances under 6 Angstrom with distance_ionic function from verif_distance module. 
    If there is not any atom of interest, an exception is raised and the interactions list is empty
    If there is not any distance according to cut-off distances, the list is also empty
    '''
    liste_ionic=recup_ionics(dico, atomes, dico_atoms)
    try :
        liste_of_ionics=distance_ionic(liste_ionic, dico_dist, dico, "IONIC", intra_inter_choice)
        liste_of_interactions=(liste_of_ionics)
    except :
        liste_of_interactions=([])      
        print("There is not any ionic bonds")
    return liste_of_interactions
'''hbond module
hbond module calculates hydrogen bonds in a protein calling hbond_call function
'''

#!/usr/bin/python3
 
from recup import *
from verif_distance_diS import *


def hbond_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice):
    '''
    hbond_call function calculates hydrogen bonds in a protein. It takes in input a list of Atome instances which constitute the protein. This list have to be created by the pdb_reader module.
    Usefull information are indicated in :
        - dico is a dictionary containing residues of interest that may participate to a specified interaction type
        - dico_atoms is a dictionary containing atoms of interest that may participate to a specified interaction type
        - dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    The module's function first gets atoms of interest with recup_accept_donr function of recup module
    Then,  it tries to calculate distances and chooses distances under 3.5 Angstrom or 4 Angstrom if a S is involved with distance_hbond function from verif_distance module.
    If there is not any atom of interest, an exception is raised and the interactions list is empty. 
    If there is not any distance according to cut-off distances, the list is also empty
    '''


    liste_tot=recup_accept_donr(dico, atomes, dico_atoms)
    liste_accept=liste_tot[0]
    liste_donr=liste_tot[1]

    try :
        liste_of_hbond=distance_hbond(liste_accept, liste_donr, dico_dist, intra_inter_choice)
        liste_of_hbond_main_main=liste_of_hbond[0]
        liste_of_hbond_main_side=liste_of_hbond[1]
        liste_of_hbond_side_side=liste_of_hbond[2]

    except :
        liste_of_hbond_main_main=[]
        liste_of_hbond_main_side=[]
        liste_of_hbond_side_side=[]     
        print("There is not any hydrogen bond")

    return (liste_of_hbond_main_main, liste_of_hbond_main_side, liste_of_hbond_side_side)


'''hydrophob module
ionic module calculates hydrophobic interactions in a protein calling hydrophob_call function
'''
#!/usr/bin/python3

from recup import *
from verif_distance_diS import *

def hydrophob_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice):
    '''
    hydroph_call function calculates hydrophobic interactions in a protein. It takes in input a list of Atome instances which constitute the protein. This list have to be created by the pdb_reader module.
    Usefull information are indicated in dico and dico_dist :
        - dico is a dictionary containing residues of interest that may participate to a specified interaction type
        - dico_atoms is a dictionary containing atoms of interest that may participate to a specified interaction type
        - dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    The module's function first gets atoms of interest with recup_hydroph function of recup module.
    Then, it tries to calculate distances and chooses distances under 5 Angstrom with distance_file function from verif_distance module. 
    If there is not any atom of interest, an exception is raised and the interactions list is empty
    If there is not any distance according to cut-off distances, the list is also empty
    '''
    liste_hydroph=recup_hydroph(dico, atomes, dico_atoms)
    try :
        liste_of_Hydroph_inter=distance_file(liste_hydroph, dico_dist, "HYDROPH", intra_inter_choice)
        liste_of_interactions=(liste_of_Hydroph_inter)
    except:
        liste_of_interactions=([])      
        print("There is not any hydrophob bonds")
    return liste_of_interactions
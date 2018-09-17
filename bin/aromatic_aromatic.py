'''aromatic_aromatic module
aromatic_aromatic module calculates aromatic/aromatic interactions in a protein calling aromatic_aromatic_call function
'''
#!/usr/bin/python3

from recup import *
from verif_distance_diS import *


def aromatic_aromatic_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice) :
    '''
    aromatic_aromatic_call function calculates aromatic/aromatic interactions in a protein. It takes in input a list of Atome instances which constitute the protein. This list have to be created by the pdb_reader module.
    Usefull information are indicated in :
        - dico is a dictionary containing residues of interest that may participate to a specified interaction type
        - dico_atoms is a dictionary containing atoms of interest that may participate to a specified interaction type
        - dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    The module's function first gets atoms of interest with recup_phenyl function of recup module.
    Then, it calculates the phenyl ring centroid, i.e. the mean position of atoms contained in the phenyl ring of a residue using centroid function from verif_distance module.
    Finaly, it tries to calculate distances and chooses distances over 4.5 Angstrom and under 7 Angstrom with distance_file function from verif_distance module.
    If there is not any atom of interest, an exception is raised and the interactions list is empty
    If there is not any distance according to cut-off distances, the list is also empty
    '''
    liste_phenyl=recup_phenyl(dico, atomes, dico_atoms)
    try :
        liste_centroid_phenyl=centroid(liste_phenyl, "phenyl") # Here, centroids are the "center" of the phenyl ring, the mean of each phenyl ring coordinates
        liste_of_phenyl_inter=distance_file(liste_centroid_phenyl, dico_dist, "PHENYL", intra_inter_choice)
        liste_of_interactions=(liste_of_phenyl_inter)
    except :
        liste_of_interactions=([])
        print("There is not any aromatic/aromatic bonds")
    return liste_of_interactions
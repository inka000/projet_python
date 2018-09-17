'''aromatic_s module
aromatic_s module calculates aromatic/S interactions in a protein calling aromatic_s_call function
'''

#!/usr/bin/python3
 
from recup import *
from verif_distance_diS import *

def aromatic_s_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice) :
    '''
    aromatic_s_call function calculates aromatic/S interactions in a protein. It takes in input a list of Atome instances which constitute the protein. This list have to be created by the pdb_reader module.
    Usefull information are indicated in :
        - dico is a dictionary containing residues of interest that may participate to a specified interaction type
        - dico_atoms is a dictionary containing atoms of interest that may participate to a specified interaction type
        - dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    The module's function first gets atoms of interest with recup_phenyl and recup_ring_S functions of recup module
    Then, it calculates the phenyl ring centroid, i.e. the mean position of atoms contained in the phenyl ring, and of the 5 members ring of TRP using centroid function from verif_distance module.
    5 members ring of TRP are considered separately
    Finaly, it tries to calculate distances and chooses distances under 5.3 Angstrom with distance_arom_s function from verif_distance module.
    If there is not any atom of interest, an exception is raised and the interactions list is empty. 
    If there is not any distance according to cut-off distances, the list is also empty
    '''

    liste_phenyl=recup_phenyl(dico, atomes, dico_atoms)
    liste_TRP_5members=recup_ring_S(dico, atomes, dico_atoms)[0]
    liste_S_MET_CYS=recup_ring_S(dico, atomes, dico_atoms)[1]

    try :
        liste_centroid_phenyl=centroid(liste_phenyl, "phenyl")
    except :
        print("There is not any phenyl ring")
        liste_centroid_phenyl=[]

    
    try :
        liste_centroid_TRP_5members=centroid(liste_TRP_5members, "5members_ring_TRP")
        liste_centroid_all=liste_centroid_phenyl+liste_centroid_TRP_5members

    except :
        liste_centroid_all=liste_centroid_phenyl
        print("There is not any TRP")
    
    try :
        liste_of_arom_s=distance_arom_s(liste_centroid_all, liste_S_MET_CYS, "AROM_S", dico_dist, intra_inter_choice)
        liste_of_interactions=(liste_of_arom_s)
    except :
        liste_of_interactions=([])
        print("There is not any aromatic/S bonds")

    return liste_of_interactions
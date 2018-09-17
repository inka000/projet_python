'''recup module
recup module allows to get specifi Atome instances from a list of Atomes instances depending on the type of interaction
'''

#!/usr/bin/python3



def recup_S(dico,liste_atomes):
    '''
    Gets the SG atoms of Cysteines residues in a list of Atome instances. 
    dico is the dictionary containing residue types corresponding to each interaction type
    Returns a list of Atome instances
    '''
    liste_S=[]
    for at in liste_atomes :
        if at.residu_type in dico["SSBOND"] and at.atome_name=="SG" :
            liste_S.append(at)
    return liste_S




def recup_hydroph(dico, liste_atomes, dico_atoms):
    '''
    Gets the atoms of hydrophobs residus in a list of Atome instances. 
    dico is the dictionary containing residue types corresponding to each interaction type
    dico_atoms is the dictionary containing atoms of interest for each interaction type
    Returns a list of Atome instances
    '''
    liste_hydroph=[]
    for at in liste_atomes :
        if at.residu_type in dico["HYDROPH"] and (at.atome_name in dico_atoms["HYDROPH"]) :
            liste_hydroph.append(at)
    return liste_hydroph



def recup_ionics(dico, liste_atomes, dico_atoms):
    '''
    Gets the atoms of ionic residus in a list of Atome instances. 
    dico is the dictionary containing residue types corresponding to each interaction type
    dico_atoms is the dictionary containing atoms of interest for each interaction type
    Returns a list of Atome instances
    ''' 
    liste_ionic=[]
    for at in liste_atomes:
        if (at.residu_type in dico["PLUS"] or at.residu_type in dico["MINUS"]) and (at.atome_name in dico_atoms["IONIC"]) :
            liste_ionic.append(at)
    return liste_ionic



def recup_phenyl(dico, liste_atomes, dico_atoms):
    '''
    Gets the atoms of aromatic residus in a list of Atome instances. 
    dico is the dictionary containing residue types corresponding to each interaction type
    dico_atoms is the dictionary containing atoms of interest for each interaction type
    This function only takes phenyl rings with 6 atoms, for the tryptophan it does not get the 5 members ring
    Returns a list of Atome instances
    ''' 
    liste_arom=[]
    for at in liste_atomes:
        if(at.residu_type in dico["PHENYL"] and at.atome_name in dico_atoms["PHENYL"]):
            if(at.residu_type=="TRP" and (at.atome_name == "CD1" or at.atome_name == "CG")): # We only want the phenyl ring
                continue
            else :
                liste_arom.append(at)
    return liste_arom




def recup_ring_S(dico, liste_atomes, dico_atoms):
    '''
    Gets the atoms of aromatic residus in a list of Atome instances. 
    dico is the dictionary containing residue types corresponding to each interaction type
    dico_atoms is the dictionary containing atoms of interest for each interaction type
    This function gets also the 5 members ring in tryptophan
    Returns two lists of Atome instances
    ''' 
    liste_5members=[]
    liste_S_MET_CYS=[]
    for at in liste_atomes :
        if (at.residu_type in dico["AROM_S"] and at.atome_name in dico_atoms["AROM_S"]) :
            liste_S_MET_CYS.append(at)
        elif (at.residu_type == "TRP" and at.atome_name in dico_atoms["AROM_TRP"]) :
            liste_5members.append(at)
    return (liste_5members, liste_S_MET_CYS)




def recup_pi(dico, liste_atomes, dico_atoms):
    '''
    Gets the atoms of lysin and arginin residus in a list of Atome instances. 
    dico is the dictionary containing residue types corresponding to each interaction type
    dico_atoms is the dictionary containing atoms of interest for each interaction type
    This function gets specific atoms that will be used to define the mean position of the cation
    Returns a list of Atome instances
    ''' 
    liste_pi=[]
    for at in liste_atomes:
        if (at.residu_type in dico["PI"] and at.atome_name in dico_atoms["PI"]):
            liste_pi.append(at)
    return liste_pi


def recup_accept_donr(dico, liste_atomes, dico_atoms):
    '''
    Gets atoms possibly involved in hydrogen bonds in a list of Atome instances. Creates two lists of Atome instances, the first contains acceptor atoms, the other contains donor atoms.
    dico is the dictionary containing residue types corresponding to each interaction type
    dico_atoms is the dictionary containing atoms of interest for each interaction type
    This function gets specific atoms that will be used to define the mean position of the cation
    Returns two lists of Atome instances
    '''
    liste_accept=[]
    liste_donr=[]
    for at in liste_atomes :
        if (at.atome_name in ["O","OXT"]): #First, the function gets all the main chain donors
            liste_accept.append(at)
        elif (at.atome_name == "N") : # Or the main chain acceptors
            liste_donr.append(at)

        if ((at.residu_type in dico["ACCEPT"] and at.residu_type in dico["DONR"]) and (at.atome_name in dico_atoms["ACCEPT"] and at.atome_name in dico_atoms["DONR"])): 
            liste_accept.append(at)
            liste_donr.append(at)
        elif (at.residu_type in dico["DONR"] and at.atome_name in dico_atoms["DONR"]) :
            liste_donr.append(at)
        elif (at.residu_type in dico["ACCEPT"] and at.atome_name in dico_atoms["ACCEPT"]) :
            liste_accept.append(at)
    return (liste_accept, liste_donr)


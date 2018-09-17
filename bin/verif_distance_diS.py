'''verif_distance module
verif_distance module allows to :
    - verifies distances between atoms or centroid in specific interaction types
    - calculate angle between two phenyl ring or between a ring and a S atom
    - calculate centroid or mean position of a phenyl ring, a 5 members ring (TRP) or an atoms group.
This module is used to caculate interactions in a protein
'''


#!/usr/bin/python3


from distance_calculation import *
from recup import *
from classes import *
from math import *



def distance_SH(liste_of_interactions, atom1, atom2, dico_dist, type_inter):
    '''
    Verifies the distance between two Atome instances provided (atom1 and atom2), based on the interaction type specified
    dico_dist is a dictionary containing distance cut-off for each interaction type
    type_inter is the interaction type (string)
    liste_of_interactions is the list of Interaction isntances the user wants to add an new interaction
    Returns the list of Interaction instances found
    '''
    #First, the function determines whether the two atoms are too far or not. 
    dist_x=abs(atom1.xpos-atom2.xpos)
    dist_y=abs(atom1.ypos-atom2.ypos)
    dist_z=abs(atom1.zpos-atom2.zpos)
    #If the two atoms are too far in one of the directions, then they cannot be nearer
    if (dist_x<=dico_dist[type_inter] or dist_y<=dico_dist[type_inter] or dist_z<=dico_dist[type_inter]) : 
        dist=distance(atom1, atom2)
        if dist>0 and dist<=dico_dist[type_inter] : #The calculated distance is verified to answer the constraints
            if (type_inter=="SHBOND" and atom1.atome_name=="SG" and atom2.atome_name=="SG" and dist<dico_dist["SSBOND"]): #Those are disulfide bridges
                pass
            else :
                liste_of_interactions.append(Interaction(atom1, atom2, dist, type_inter))
    return liste_of_interactions





def angle(atom1, atom2, type_inter):
    '''
    Calculates the angle between two phenyl ring or between a phenyl ring and a S atom.
    atom1 and atom2 are Atome insances. If the interaction is aroamtic/S, atom1 must be the centroid and atom2 the S atom
    type_inter is the interaction type (string)
    Calculates the normal vector of a plan with the cross product
    The first Atome instance must be a centroid of an aromatic ring with, instead of the atome_num, two Atome instances included in this ring
    The angle will be between the normal vector and the vector between the centroid of the other ring and an atom of this ring (aromatic/aromatic), or with the centroid of the first ring (aromatic/S)
    Returns the angle in degrees with the functions acos and degrees of math module
    '''
    #First, we will calculate the normal vector of the phenyl ring
    vec1=[atom1.xpos-atom1.atome_num[0].xpos, atom1.ypos-atom1.atome_num[0].ypos, atom1.zpos-atom1.atome_num[0].zpos] # vector between the first centroid and an atom of the phenyl ring associated with
    vec2=[atom1.xpos-atom1.atome_num[1].xpos, atom1.ypos-atom1.atome_num[1].ypos, atom1.zpos-atom1.atome_num[1].zpos] #vector between the first centroid and the second atom of the phenyl ring
    cross=[(vec1[1]*vec2[2])-(vec1[2]*vec2[1]),(vec1[2]*vec2[0])-(vec1[0]*vec2[2]),(vec1[0]*vec2[1])-(vec1[1]*vec2[0])] #Vector of the cross product (normal to the plan)
    # The normal vector calculated is not oriented, thus the angle calculated can be 90 too big

    if type_inter =="PHENYL": #the vector between the centroid of the other phenyl ring and an atom of this ring is calcualted
        vec3=[atom2.xpos-atom2.atome_num[0].xpos, atom2.ypos-atom2.atome_num[0].ypos, atom2.zpos-atom2.atome_num[0].zpos]
    else :
        vec3=[atom1.xpos-atom2.xpos, atom1.ypos-atom2.ypos, atom1.zpos-atom2.zpos] # the vector between the S atom and the centroid of the phenyl ring
    x=((cross[0]*vec3[0])+(cross[1]*vec3[1])+(cross[2]*vec3[2]))/(sqrt(((cross[0]**2)+(cross[1]**2)+(cross[2]**2))*((vec3[0]**2)+(vec3[1]**2)+(vec3[2]**2))))
    ang=acos(x)
    return degrees(ang)




def distance_phenyl(liste_of_interactions, atom1, atom2, dico_dist, type_inter):
    '''
    Verifies the distance between two phenyl ring centroids (atom1 and atom2 which are Atome instances)
    liste_of_interactions is the list of Interaction isntances the user wants to add an new interaction
    type_inter is the interaction type (string)
    dico_dist is a dictionary containing distance cut-off for each interaction type
    Returns a list of Interaction instances
    '''
    dist_x=abs(atom1.xpos-atom2.xpos)
    dist_y=abs(atom1.ypos-atom2.ypos)
    dist_z=abs(atom1.zpos-atom2.zpos)
    verif_x=(dist_x<dico_dist[type_inter][1]) #If the atom is to far, it cannot be brought closer, however, the opposit is possible
    verif_y=(dist_y<dico_dist[type_inter][1])
    verif_z=(dist_z<dico_dist[type_inter][1])

    if (verif_x and verif_y and verif_z) :
        dist=distance(atom1, atom2)
        ang=angle(atom1, atom2, type_inter) #The angle is an information
        if (dist > dico_dist[type_inter][0] and dist < dico_dist[type_inter][1]) : # It verifies that the distance answers the constraints. 
            liste_of_interactions.append(Interaction(atom1, atom2, dist, type_inter, ang ))

    return liste_of_interactions




def distance_phenyl_S(liste_of_interactions, atom1, atom2, dico_dist, type_inter):
    '''
    Verifies the distance between one ring centroid (atom1) and one S atom (atom2) which are Atome instances
    liste_of_interactions is the list of Interaction isntances the user wants to add an new interaction
    type_inter is the interaction type (string)
    dico_dist is a dictionary containing distance cut-off for each interaction type
    Returns a list of Interaction instances
    '''
    dist_x=abs(atom1.xpos-atom2.xpos)
    dist_y=abs(atom1.ypos-atom2.ypos)
    dist_z=abs(atom1.zpos-atom2.zpos)
    verif_x=(dist_x<dico_dist[type_inter]) #If the atom is to far, it cannot be brought closer, however, the opposit is possible
    verif_y=(dist_y<dico_dist[type_inter])
    verif_z=(dist_z<dico_dist[type_inter])

    if (verif_x and verif_y and verif_z) :
        dist=distance(atom1, atom2)
        ang=angle(atom1, atom2, type_inter) #The angle is an information
        if (dist > 0 and dist < dico_dist[type_inter]): #the function has to verify that the distance is not too small
            liste_of_interactions.append(Interaction(atom1, atom2, dist, type_inter, ang))

    return liste_of_interactions




def distance_file(liste_atom, dico_dist, type_inter, intra_inter_choice) :
    '''
    Creates a list of interactions
    liste_atom contains the Atome instances of interest
    Interaction type must be specified with type_inter (string), the distance calculation function used depends on the interaction
    dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    Returns a list of Interaction instances
    '''
    liste_of_interactions=[]
    actual_atom = 0
    while actual_atom<(len(liste_atom)-1) : #The actual atom is used with all next atoms, not without itself and not with previous atoms since it already has been done
        for i in range((actual_atom+1), len(liste_atom)) :
            if(liste_atom[actual_atom].residu_num == liste_atom[i].residu_num and liste_atom[actual_atom].chain == liste_atom[i].chain): #the atoms come from different residues
                pass
            else :
                if (intra_inter_choice==0) :
                    if liste_atom[actual_atom].chain == liste_atom[i].chain : #The atoms have to be on the same chain
                        if type_inter!="PHENYL" :
                            liste_of_interactions=distance_SH(liste_of_interactions, liste_atom[actual_atom], liste_atom[i], dico_dist, type_inter) 
                        else :
                            liste_of_interactions=distance_phenyl(liste_of_interactions, liste_atom[actual_atom], liste_atom[i], dico_dist, type_inter )
                else :
                    if liste_atom[actual_atom].chain != liste_atom[i].chain : #The atoms cannot be on the same chain
                        if type_inter!="PHENYL" :
                            liste_of_interactions=distance_SH(liste_of_interactions, liste_atom[actual_atom], liste_atom[i], dico_dist, type_inter) 
                        else :
                            liste_of_interactions=distance_phenyl(liste_of_interactions, liste_atom[actual_atom], liste_atom[i], dico_dist, type_inter )
        actual_atom +=1 
    return liste_of_interactions


def distance_ionic(liste_atom, dico_dist, dico, type_inter, intra_inter_choice):
    '''
    Creates a list of ionic interactions 
    liste_atom contains the Atome instances of interest
    Interaction type must be specified with type_inter (string), the distance calculation function used depends on the interaction
    dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    Returns a list of Interaction instances
    '''
    liste_of_interactions=[]
    actual_atom=0
    while actual_atom<(len(liste_atom)-1) :
        for i in range((actual_atom+1), len(liste_atom)) : #The actual atom is used with all next atoms, not without itself and not with previous atoms since it already has been done
            if (intra_inter_choice==0) :
                if liste_atom[actual_atom].chain == liste_atom[i].chain : #The atoms have to be on the same chain
                    #atoms must have opposite charge
                    if(liste_atom[actual_atom].residu_type in dico["PLUS"] and  liste_atom[i].residu_type in dico["MINUS"]) or (liste_atom[actual_atom].residu_type in dico["MINUS"] and liste_atom[i].residu_type in dico["PLUS"]):
                        liste_of_interactions=(distance_SH(liste_of_interactions, liste_atom[actual_atom], liste_atom[i], dico_dist, type_inter ))
            else :
                if liste_atom[actual_atom].chain != liste_atom[i].chain : #The atoms cannot be on the same chain
                    #atoms must have opposite charge
                    if(liste_atom[actual_atom].residu_type in dico["PLUS"] and  liste_atom[i].residu_type in dico["MINUS"]) or (liste_atom[actual_atom].residu_type in dico["MINUS"] and liste_atom[i].residu_type in dico["PLUS"]):
                        liste_of_interactions=(distance_SH(liste_of_interactions, liste_atom[actual_atom], liste_atom[i], dico_dist, type_inter ))
        actual_atom +=1
    return liste_of_interactions


def centroid(liste_atom, name):
    '''
    Calculates the centroid of an aromatic ring or the mean position of an atoms group
    Takes a liste of Atome instances
    The name specified will be the name instead of atome.name attribute
    Return a list of centroid (Atome instances)
    '''
    s_x=liste_atom[0].xpos; s_y=liste_atom[0].ypos; s_z=liste_atom[0].zpos #First, the function get the first atom of the first ring coordinates
    cpt_atoms=1 #The atom counter starts at 1
    liste_centroid=[]
    for i in range(1, len(liste_atom)) :
        if(liste_atom[i].residu_num == liste_atom[i-1].residu_num) and (liste_atom[i].chain == liste_atom[i-1].chain) : #If this is true, the atom is in the same residue and chain that the previous atom
            cpt_atoms+=1 #It is then an other atom of the same phenyl ring
            s_x+=liste_atom[i].xpos
            s_y+=liste_atom[i].ypos
            s_z+=liste_atom[i].zpos
        else : #The phenyl ring is complete
            # Instead of atom number, we keep two atoms of the ring. This will be used for angle calculation.
            liste_centroid.append(Atome(name, liste_atom[i-1].chain, s_x/cpt_atoms, s_y/cpt_atoms, s_z/cpt_atoms, liste_atom[i-1].residu_type, liste_atom[i-1].residu_num, [liste_atom[i-1], liste_atom[i-2]]))
            # i is on the next atom which is not included in the ring but in the next one
            cpt_atoms=1
            s_x=liste_atom[i].xpos # positions are reset to the actual i atom
            s_y=liste_atom[i].ypos
            s_z=liste_atom[i].zpos
    liste_centroid.append(Atome(name, liste_atom[i-1].chain, s_x/cpt_atoms, s_y/cpt_atoms, s_z/cpt_atoms, liste_atom[i-1].residu_type, liste_atom[i-1].residu_num, [liste_atom[i-1], liste_atom[i-2]]))

    return liste_centroid



def distance_arom_s(liste_centroid, liste_S, type_inter, dico_dist, intra_inter_choice):
    '''
    Creates a list of aromatic/S interactions using S atoms from methionine and cysteine residue and aromatic ring centroids
    Interaction type must be specified, the distance calculation function used depends on the interaction
    dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    Returns a list of Interaction instances
    '''
    liste_of_interactions=[]
    for i in range(len(liste_S)):
        for j in range(len(liste_centroid)): #Each S atom will be tested with each phenyl ring centroid of the same chain
            if (intra_inter_choice==0) :
                if liste_S[i].chain == liste_centroid[j].chain :  #The atoms have to be on the same chain
                    liste_of_interactions=(distance_phenyl_S(liste_of_interactions, liste_centroid[j], liste_S[i], dico_dist, type_inter))
            else :
                if liste_S[i].chain != liste_centroid[j].chain : #The atoms cannot be on the same chain
                    liste_of_interactions=(distance_phenyl_S(liste_of_interactions, liste_centroid[j], liste_S[i], dico_dist, type_inter))
            
    return liste_of_interactions



def distance_hbond(liste_accept, liste_donr, dico_dist, intra_inter_choice):
    '''
    Creates a list of hydrogen interactions
    Depending on the atom involved, the interaction can be main chain/main chain, main chain/side chain or side chain/side chain.
    liste_accept is a list of Atomes instances considered as acceptors
    liste_donr is a list of Atomes instances considered as donors
    dico_dist is a dictionary containing distance cut-off for each interaction type
    intra_inter_choice is the choice between intra-protein and protein-protein intraction (integer): 0 for intra; 1 for protein-protein
    Returns three list of Interaction instances
    '''
    liste_of_interactions_main_main=[]
    liste_of_interactions_main_side=[]
    liste_of_interactions_side_side=[]
    for i in range(len(liste_donr)):
        for j in range(len(liste_accept)):
            if (intra_inter_choice==0) :
                if (liste_accept[j].chain == liste_donr[i].chain) :   #The atoms have to be on the same chain
                    if (liste_accept[j].residu_num<(liste_donr[i].residu_num-1) or liste_accept[j].residu_num>(liste_donr[i].residu_num+1)):
                        if liste_accept[j].atome_name in ["SD","SG"] or liste_donr[i].atome_name in ["SD","SG"] :
                            type_inter="SHBOND"
                        else :
                            type_inter="HBOND"
                        if (liste_accept[j].atome_name in ["O","OXT"] and liste_donr[i].atome_name == "N"):
                            liste_of_interactions_main_main=distance_SH(liste_of_interactions_main_main, liste_donr[i], liste_accept[j], dico_dist, type_inter)
                        elif (liste_accept[j].atome_name in ["O","OXT"] or liste_donr[i].atome_name == "N"):
                            liste_of_interactions_main_side=distance_SH(liste_of_interactions_main_side, liste_donr[i], liste_accept[j], dico_dist, type_inter)
                        else :
                            liste_of_interactions_side_side=distance_SH(liste_of_interactions_side_side, liste_donr[i], liste_accept[j], dico_dist, type_inter)
            else :
                if (liste_accept[j].chain != liste_donr[i].chain) : #The atoms cannot be on the same chain
                    if (liste_accept[j].residu_num<(liste_donr[i].residu_num-1) or liste_accept[j].residu_num>(liste_donr[i].residu_num+1)):
                        if liste_accept[j].atome_name in ["SD","SG"] or liste_donr[i].atome_name in ["SD","SG"] :
                            type_inter="SHBOND"
                        else :
                            type_inter="HBOND"
                        if (liste_accept[j].atome_name in ["O","OXT"] and liste_donr[i].atome_name == "N"):
                            liste_of_interactions_main_main=distance_SH(liste_of_interactions_main_main, liste_donr[i], liste_accept[j], dico_dist, type_inter)
                        elif (liste_accept[j].atome_name in ["O","OXT"] or liste_donr[i].atome_name == "N"):
                            liste_of_interactions_main_side=distance_SH(liste_of_interactions_main_side, liste_donr[i], liste_accept[j], dico_dist, type_inter)
                        else :
                            liste_of_interactions_side_side=distance_SH(liste_of_interactions_side_side, liste_donr[i], liste_accept[j], dico_dist, type_inter)
    return (liste_of_interactions_main_main, liste_of_interactions_main_side, liste_of_interactions_side_side)
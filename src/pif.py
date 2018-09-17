'''
Source program of the protein calculation interaction software.
Allows to calculate different types of interaction in a protein in pdb format :
    - Disulfide bridges
    - ionic interactions
    - hydrophobic interactions
    - aromatic/aromatic interactions
    - aromatic/S interactions
    - aromatic/cation-pi interactions
    - hydrogen bonds
The user will be asked to choose the type(s) of interactions and whether he/she wants intra-protein or protein-protein interactions. 
For protein-protein interaction, the pdb file provided must have multiple chains.
To use this program, the user must type :
    python3 pif.py file.pdb
Exemple :
    >>python3 pif.py 1zni.pdb
    >>1
    SS bridges of the protein 1zni will be calculated and the file 1zni.txt will be created with tabulation between column :

    SS bridges

    Position    Residue Chain   Position    Residue Chain   Distance(A) Angle
    6   CYS A   11  CYS A   2.031   999.900
    7   CYS A   7   CYS B   2.023   999.900
    20  CYS A   19  CYS B   2.067   999.900
    6   CYS C   11  CYS C   1.990   999.900
    7   CYS C   7   CYS D   2.022   999.900
    20  CYS C   19  CYS D   2.067   999.900
    999.9 angle means that the angle was not calculated


'''

#!/usr/bin/python3

import sys
sys.path.insert(0,'../bin') #the path of the modules is set to the bin/ directory
from pdb_reader import *
from write_txt import *
from classes import *
from dicos import *
from ssbridges import *
from ionic import *
from hydrophob import *
from aromatic_aromatic import *
from aromatic_s import *
from pi_cation_aromatic import *
from hbond import *

def main():
    '''
    Main function of the program. It creates lists of atoms, verifies that the file provided does exist in the data directory
    The user is asked which interaction type(s) he/she wants to calculate. If the user does not provide an answer, the program stops with an error
    All interactions are tested but when they do not exist the program does not stop, it continues with the next interaction type
    
    The file studied must be an argument and must exist, if not, a message indicates that it does not. 
    The program creates a txt file with all informations about interaction types specified by the user. The file is named according to the pdb name
    '''
    try :
       file=sys.argv[1]
    except IndexError :
       print("There is not any file provided, please provide a file name\n")
    else : 
        try : 
            atomes=readPDB(file)
        except IOError :
            print("The file does not exist in the directory, please provide an existing file\n")  
        else :      
            response = list(input("Choose the type of interactions you want to find typing the corresponding numbers without space \n 1 : SS bridges\n 2: ionic interactions\n 3: hydrophobes interactions\n 4: aromatic/aromatic interactions\n 5: aromatic/S interactions\n 6: aromatic/pi interactions\n 7: hydrogen interactions\n"))
            assert(len(response)!=0 and len(response)<8) #If the assertion is false, the program stops with an error message
            rep=[]
            for num in response :
                assert(int(num)<8 and int(num)>0)
                rep.append(int(num))

            intra_inter_choice = int(str(input("Choose if the interactions have to be intra-protein [0] or protein-protein [1] :")))
            assert(intra_inter_choice==0 or intra_inter_choice==1) #If the assertion is false, the program stops with an error message
            
            liste_of_interactions={} #This is the dictionary of all interactions calculated and found


            if 1 in rep :
                liste_of_interactions["SS bridges"]=ssbridge_call(dico, atomes, dico_dist, intra_inter_choice)

            if 2 in rep :
                liste_of_interactions["ionic interactions"]=ionic_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice)
         
            if 3 in rep :
                liste_of_interactions["hydrophob interactions"]=hydrophob_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice)
                
            if 4 in rep :
                liste_of_interactions["aromatic/aromatic interactions"]=aromatic_aromatic_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice)


            if 5 in rep :
                liste_of_interactions["aromatic/S interactions"]=aromatic_s_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice)
                
            if 6 in rep :
                liste_of_interactions["pi/aromatic interactions"]=pi_cation_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice)

            if 7 in rep :
                liste_tot=hbond_call(dico, atomes, dico_dist, dico_atoms, intra_inter_choice)
                liste_of_interactions["HBOND main_main"]=liste_tot[0]
                liste_of_interactions["HBOND main_side"]=liste_tot[1]
                liste_of_interactions["HBOND side_side"]=liste_tot[2]

            write_file_resume(liste_of_interactions, file.split(".")[0])    
        
    return 0
    

if __name__=="__main__":
    main()

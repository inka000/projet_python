'''write_txt module
write_txt module allows to write two txt files resuming interactions information with wrtie_file_resume function
'''

#!/usr/bin/python3

from classes import *


def write_file_resume(liste, name):
    '''
    This function write two txt files. The user must provide a dictionnary of Interactions class instances. The keys are the names of interactions. 
    The user must provide the name of the created file
    The first file is user friendly with informations that the user can read easily.
    The other file resume all residues involved in interactions and can be read for example with R software
    Each column is separated with the next one with a tabulation
    '''
    f2=open("../results/"+name+".txt", "w")
    f=open("../results/"+"read"+name+".txt", "w")
    for inter in liste :
        f.write(str(inter))
        f.write("\n\n")
        if (inter in ["SS bridges","aromatic/aromatic interactions","aromatic/S interactions","pi/aromatic interactions"]) :
            f.write("Position\tResidue\tChain\tPosition\tResidue\tChain\tDistance(A)\tAngle")
            f.write("\n")
            for i in range(len(liste[inter])):
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}".format(liste[inter][i].atom1.residu_num, liste[inter][i].atom1.residu_type, liste[inter][i].atom1.chain, liste[inter][i].atom2.residu_num, liste[inter][i].atom2.residu_type, liste[inter][i].atom2.chain, liste[inter][i].dist, liste[inter][i].angle ))
                f.write("\n")
                f2.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}".format(liste[inter][i].atom1.residu_num, liste[inter][i].atom1.residu_type, liste[inter][i].atom1.chain, liste[inter][i].atom2.residu_num, liste[inter][i].atom2.residu_type, liste[inter][i].atom2.chain, liste[inter][i].dist, liste[inter][i].angle ))
                f2.write("\n")
        else : 
            f.write("Position\tAtom\tResidue\tChain\tPosition\tAtom\tResidue\tChain\tDistance(A)\tAngle")
            f.write("\n")
            for i in range(len(liste[inter])):
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}".format(liste[inter][i].atom1.residu_num, liste[inter][i].atom1.atome_name, liste[inter][i].atom1.residu_type, liste[inter][i].atom1.chain, liste[inter][i].atom2.residu_num, liste[inter][i].atom2.atome_name, liste[inter][i].atom2.residu_type, liste[inter][i].atom2.chain, liste[inter][i].dist, liste[inter][i].angle ))
                f.write("\n")
                f2.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}".format(liste[inter][i].atom1.residu_num, liste[inter][i].atom1.residu_type, liste[inter][i].atom1.chain, liste[inter][i].atom2.residu_num, liste[inter][i].atom2.residu_type, liste[inter][i].atom2.chain, liste[inter][i].dist, liste[inter][i].angle ))
                f2.write("\n")
        f.write("999.9 angle means that the angle was not calculated\n\n")
    f.close()

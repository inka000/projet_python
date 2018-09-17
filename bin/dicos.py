'''dicos module
This module contains :
    - dico : a dictionary of residues of interest depending on the interaction type
    - dico_dist : a dictionary of distance cut-offs depending on the interaction type
    - dico_atoms : a dictionary of atoms of interest depending on the interaction type
'''


#!/usr/bin/python3


dico={}
dico["SSBOND"]=["CYS"]
dico["HYDROPH"]=["ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","TYR"]
dico["PLUS"]=["ARG","LYS","HIS"] # positively charged residues
dico["MINUS"]=["GLU", "ASP"] # negatively charged residues
dico["PHENYL"]=["PHE","TYR","TRP"]
dico["AROM_S"]=["CYS","MET"]
dico["PI"]=["LYS","ARG"]
dico["ACCEPT"]=["MET","ASN","GLN","CYS","ASP","GLU","HIS","TYR","SER","THR"]
dico["DONR"]=["ARG","LYS","TRP","ASN","GLN","CYS","ASP","GLU","HIS","TYR","SER","THR"]

dico_dist={}
dico_dist["SSBOND"]=2.2
dico_dist["HYDROPH"]=5.0
dico_dist["IONIC"]=6.0
dico_dist["PHENYL"]=[4.5, 7.0]
dico_dist["AROM_S"]=5.3
dico_dist["PI"]=6.0
dico_dist["HBOND"]=3.5
dico_dist["SHBOND"]=4.0 #cut-off of hydrogen bond when a S is involved

dico_atoms={}
dico_atoms["HYDROPH"]=["CB","CG","CB1","CG2","CD","CD1","CD2","CE","CE1","CE2","CZ","CZ2","CZ3","CH2","NE1","OH"]
dico_atoms["IONIC"]=["OE1","OE2","OD1","OD2","NE","NH1","NH2","NZ","NE2","ND1"]
dico_atoms["PHENYL"]=["CG","CD1","CD2","CE1","CE2","CZ","CE3","CZ2","CZ3","CH2"]
dico_atoms["AROM_TRP"]=["CG","CD1","CD2","CE2", "NE1"] # 5 members TRP ring 
dico_atoms["AROM_S"]=["SG","SD"] # S atoms involved in aromatic/S interactions
dico_atoms["PI"]=["NE","CZ","NH1","NH2","CE","NZ"]
dico_atoms["ACCEPT"]=["SD","OD1","OD2","OE1","OE2","OH","NE2","ND2","ND1","SG","OG","OG1"]
dico_atoms["DONR"]=["NH1","NH2","NZ","NE1","OD1","OD2","OE1","OE2","OH","NE2","ND2","ND1","SG","OG","OG1"]
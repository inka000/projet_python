PIF (Protein Interaction Finder) is a software that can find interactions within a pdb file 
of a protein or a complex of proteins
To use PIF the user must open a Linux terminal and type : python3 pif.py filename.pdb
Note that pif.py is in src/ directory and the pdb file must be in the data/ directory.

The user will be asked which type(s) of interactions he/she wants to find. 
The user has to give an answer typing at least one possibility. Response is the numbers corresponding to
interaction types, without space. 
Then, the user must type 0 for intra-protein interactions or 1 for inter-proteins interactions.
Both at the same time is not possible. For inter-proteins interactions, the pdb file must contains at least 
two different chains. 
PIF is based on official pdb notations and format, the pdb file has to be an official pdb file. 

It is possible to change distance cut_offs for each interaction. To do so, please open dicos.py to edit it. 
dico_dist contains the distances. 
''
dico_dist["SSBOND"]=2.2
dico_dist["HYDROPH"]=5.0
dico_dist["IONIC"]=6.0
dico_dist["PHENYL"]=[4.5, 7.0]
dico_dist["AROM_S"]=5.3
dico_dist["PI"]=6.0
dico_dist["HBOND"]=3.5
dico_dist["SHBOND"]=4.0
''
SSBOND : disulfide bridges
HYDROPH : hydrophobic interactions
IONIC : ionic interactions
PHENYL : aromatic/aromatic interactions (min and max)
AROM_S : aromatic/sulfur interactions
PI : cation-pi interactions
HBOND : hydrogen bonds without a sulfur atom involved
SHBOND : hydrogen bonds with a sulfur atom involved

Once interaction types and inter/intra specification are chosen, PIF will create two .txt files :
- readfilename.txt which contains information about interactions that the user can easily read
- filename.txt which contains all residues involved in an interaction

Interactions calculation criteria :
Disulfide bridges : distance between two cysteine sulfur atom is calculated. There is an interaction when the 
distance is lower than 2.2 A. 

Hydrophobic interactions : distance between two hydrophobic groups atoms is calculated. There is an interaction 
when the distance is lower than 5 A. Alanine, valine, leucine, isoleucine, methionine, phenylalanine, tryptophan, 
proline, tyrosine are considered as hydrophobic amino acids. 

Ionic interactions : distance between opposite charged atoms of ionic residues is calculated. When the distance 
is under 6 A there is an interaction. Asparagine and glutamate are considered as negative residues, and lysine, 
histidine and arginine are considered as positive residues. 

Hydrogen interactions : distances between pairs of potential hydrogen acceptors or donors are calculated. The 
distance cut-off depends on the atoms involved. If a sulfur is involved, the distance must be under 4 A, 
otherwise, the distance must be under 3.5 A. The selection of acceptors or donors depend on JOY software 
criteria. Hydrogen bonds can be between main chains, main chain and side chain, or side chains. For SG atoms, 
if the distance between two SG atoms is lower than the disulfide bridge cut-off, it is not a hydrogen bond but 
a disulfide bridge.

Aromatic/aromatic interactions : phenyl ring of tyrosine, phenylalanine and tryptophane are considered. First, 
centroids are determined, i.e. the ring center coordinates. Then, distances between pairs of centroids are 
calculated. For an interaction, the distance must be between 4.5 and 7 A. When distances are determined, 
the coordinates of the normal vector of the first phenyl ring are calculated with a cross product of two 
non collinear vectors collinear with the ring. The angle between this vector and a collinear vector of the 
second phenyl ring is calculated.

Aromatic/Sulfur interactions :phenyl ring of tyrosine, phenylalanine and tryptophane are considered, as well 
as the 5 members ring of tryptophane separately. Centroid are calculated like in aromatic/aromatic interactions. 
The other participant is a sulfur atom from a cysteine or a methionine. The distance between the sulfur atom 
and the ring centroid must be lower than 5.3 A for an interaction. The normal vector of the ring is determined 
like in the aromatic/aromatic interactions. The coordinates of the vector between the centroid and the sulfur 
atom are calculated and the angle between the normal vector and this last vector.

Cation-π interactions : according to R. Sathyapriya and S. Vishveshwara, a cation-π interaction involved a 
positive charge center of an arginine or a lysine and an aromatic ring. Like in aromatic/sulfur interactions, 
the 5 members ring of tryptophane is considered separately. Also, the positive charge centers are the mean 
coordinates of NE, CZ, NH1 and NH2 for arginine, CE and NZ for lysine. The distance must be under 6 A. 

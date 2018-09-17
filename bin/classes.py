'''classes module
classes modules manages Atome and Interaction classes
'''
#!/usr/bin/python3



class Atome :
    '''
    class Atome : object with atoms informations
    the name of the atom
    the chain the atom comes from
    x,y and z coordinates of the atom
    residu type with 3 letters code name
    residu number in the protein
    atom number in the protein
    This class is used for centroid calculation because it needs the same informations
    '''
    def __init__(self, atome_name, chain, xpos, ypos, zpos, residu_type, residu_num, atome_num):
        '''
        __init__ method initiates an Atome instance
        '''
        self.atome_name=atome_name
        self.chain=chain
        self.xpos=xpos
        self.ypos=ypos
        self.zpos=zpos
        self.residu_type=residu_type 
        self.residu_num=residu_num 
        self.atome_num=atome_num
    def __str__(self):
        '''
        __str__ method can be used to print an Atome instance attributes
        '''
        return "name :{} chain :{} x :{} y :{} z :{} residu_type :{} residu_num :{} atome_num :{}".format(self.atome_name, self.chain, self.xpos, self.ypos , self.zpos, self.residu_type, self.residu_num, self.atome_num)




class Interaction:
    '''
    class Interaction : instance with interaction informations
    the number of interactions in total
    the first atom in the interaction (instance of Atome class)
    the second atom in the interaction (instance of Atome class)
    the distance of the interaction
    the angle of the interaction, 0.0 by default
    ''' 
    numOfInstances = 0
    def __init__(self, atom1, atom2, dist, type_inter, angle = 999.9):
        '''
        __init__ method initiates an Interaction instance
        '''
        Interaction.numOfInstances += 1
        self.atom1=atom1
        self.atom2=atom2
        self.dist=dist
        self.type_inter=type_inter
        self.angle=angle
        




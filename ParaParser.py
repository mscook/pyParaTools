import sys
from Bio.PDB import *

class ParaParser:
    def __init__(self, stdin):
        self.data_type    = stdin[1]
        self.pdb_fn       = stdin[2]
        self.para_data_fn = stdin[3]
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self.pdb_fn)
        self.structure    = structure_in
        self.model        =self.structure[0]
        para_data_in = open(self.para_data_fn).readlines()
        self.dataset      = para_data_in
        self.parsed       = []

    def getDataType(self):
        return self.data_type

    def getPDBFn(self):
        return self.pdb_fn

    def getParaDataFn(self):
        return self.para_data_fn

    def getStructure(self):
        return self.structure

    def getDataset(self):
        return self.dataset

    def getNumModels(self):
        nmodels = 0
        for model in self.structure:
            nmodels = nmodels+1
        return nmodels

    def getParsed(self):
        return self.parsed

    def setParaDataFn(self, data_name):
        self.data_type = data_name

    def setPDBFn(self, pdb_name):
        self.pdb_fn = pdb_name

    def setParaDataFn(self, para_name):
        self.para_data_fn = para_name

    def setModel(self, num):
        self.model = self.structure[num]

    def doParse(self):
        results = []
        for i in range(0, len(self.dataset)):
            resi, atom_type, exp, tol =  self.dataset[i].split()
            for atom in self.model.get_atoms():
                cur_resi      = str(list(atom.get_parent().get_id())[1])
                cur_atom_type = atom.get_name()
                if (cur_resi.strip() == resi.strip()) and (cur_atom_type.strip()
                 == atom_type.strip()):
                    tmp = []
                    tmp.append(atom_type.strip())
                    tmp.append(resi.strip())
                    tmp.append(exp)
                    tmp.append(tol)
                    tmp.append(atom.get_coord().tolist())
                    results.append(tmp)
        self.parsed = results



class PCSParser(ParaParser):
    def __init__(self, stdin):
        self.data_type    = stdin[1]
        self.pdb_fn       = stdin[2]
        self.para_data_fn = stdin[3]
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self.pdb_fn)
        self.structure    = structure_in
        self.model        =self.structure[0]
        para_data_in = open(self.para_data_fn).readlines()
        self.dataset      = para_data_in
        self.metal_loc    = []
        self.metal_loc.append(float(sys.argv[4]))
        self.metal_loc.append(float(sys.argv[5]))
        self.metal_loc.append(float(sys.argv[6]))
        self.ax           = stdin[7]
        self.rh           = stdin[8]
        self.Ealpha       = stdin[9]
        self.Ebeta        = stdin[10]
        self.Egamma       = stdin[11]



    def getMetalLoc(self):
        return self.metal_loc

    def getMetalLocx(self):
        return self.metal_loc[1]

    def getMetalLocy(self):
        return self.metal_loc[2]

    def getMetalLocz(self):
        return self.metal_loc[2]

    def getAxial(self):
        return self.ax

    def getRhombic(self):
        return self.rh

    def getAlpha(self):
        return self.Ealpha

    def getBeta(self):
        return self.Ebeta

    def getGamma(self):
        return self.Egamma


    def setMetalLoc(self, metal_xyz):
        self.metal_loc = metal_xyz

    def setMetalLocx(self, mx):
        self.metal_loc[1] = mx

    def setMetalLocy(self, my):
        self.metal_loc[2] = my

    def setMetalLocz(self, mz):
        self.metal_loc[2] = mz

    def setAxial(self, axial):
        self.ax = axial

    def setRhombic(self, rhombic):
        self.rh = rhombic

    def setAlpha(self, alpha):
        self.Ealpha = alpha

    def setBeta(self, beta):
        self.Ebeta = beta

    def setGamma(self, gamma):
        self.Egamma = gamma



class PREParser(ParaParser):
    def __init__(self, stdin):
        self.data_type    = stdin[1]
        self.pdb_fn       = stdin[2]
        self.para_data_fn = stdin[3]
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self.pdb_fn)
        self.structure    = structure_in
        self.model        =self.structure[0]
        para_data_in = open(self.para_data_fn).readlines()
        self.dataset      = para_data_in
        self.metal_loc    = []
        self.metal_loc.append(float(sys.argv[4]))
        self.metal_loc.append(float(sys.argv[5]))
        self.metal_loc.append(float(sys.argv[6]))
        self.c            = float(stdin[7])

    def getMetalLoc(self):
        return self.metal_loc

    def getMetalLocx(self):
        return self.metal_loc[1]

    def getMetalLocy(self):
        return self.metal_loc[2]

    def getMetalLocz(self):
        return self.metal_loc[2]

    def getConstant(self):
        return self.c

    def setMetalLoc(self, metal_xyz):
        self.metal_loc = metal_xyz

    def setMetalLocx(self, mx):
        self.metal_loc[1] = mx

    def setMetalLocy(self, my):
        self.metal_loc[2] = my

    def setMetalLocz(self, mz):
        self.metal_loc[2] = mz

    def setConstant(self, constant):
        self.c = constant



class RDCParser(ParaParser):
    def __init__(self, stdin):
        self.data_type    = stdin[1]
        self.pdb_fn       = stdin[2]
        self.para_data_fn = stdin[3]
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self.pdb_fn)
        self.structure    = structure_in
        self.model        =self.structure[0]
        para_data_in = open(self.para_data_fn).readlines()
        self.dataset      = para_data_in
        self.ax           = stdin[4]
        self.rh           = stdin[5]
        self.Ealpha       = stdin[6]
        self.Ebeta        = stdin[7]
        self.Egamma       = stdin[8]

    def getAxial(self):
        return self.ax

    def getRhombic(self):
        return self.rh

    def getAlpha(self):
        return self.Ealpha

    def getBeta(self):
        return self.Ebeta

    def getGamma(self):
        return self.Egamma


    def setAxial(self, axial):
        self.ax = axial

    def setRhombic(self, rhombic):
        self.rh = rhombic

    def setAlpha(self, alpha):
        self.Ealpha = alpha

    def setBeta(self, beta):
        self.Ebeta = beta

    def setGamma(self, gamma):
        self.Egamma = gamma

    def doParse(self):
        results = []
        for i in range(0, len(self.dataset)):
            resi, atom_type, exp, tol =  self.dataset[i].split()
            at1, at2 = atom_type[1], atom_type[0]
            for atom in self.model.get_atoms():
                cur_resi      = str(list(atom.get_parent().get_id())[1])
                cur_atom_type = atom.get_name()
                c_at1, c_at2 = cur_atom_type[1], cur_atom_type[0]
                tmp1 = []
                if (cur_resi.strip() == resi.strip()) and (at1.strip()
                 == c_at1.strip()):
                    tmp.append(atom_type.strip())
                    tmp.append(resi.strip())
                    tmp.append(exp)
                    tmp.append(tol)
                    tmp.append(atom.get_coord().tolist())
                if (cur_resi.strip() == resi.strip()) and (at2.strip()
                 == c_at2.strip()):


                    results.append(tmp)
        self.parsed = results
        return atom1, atom2


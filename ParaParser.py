import sys, math
from Bio.PDB import *
from ParaData import *

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
        #This method (and the RDC one) needs a bit of an optimize.
        pDlist = []
        for i in range(0, len(self.dataset)):
            resi, atom_type, exp, tol =  self.dataset[i].split()
            for atom in self.model.get_atoms():
                cur_resi      = str(list(atom.get_parent().get_id())[1])
                cur_atom_type = atom.get_name()
                if (cur_resi.strip() == resi.strip()) and (cur_atom_type.strip()
                 == atom_type.strip()):
                    #This is possibly a hack...
                    if self.data_type.strip() == 'pcs':
                        pDlist.append(PCSData(atom_type.strip(), resi.strip(),
                        exp, tol, atom.get_coord().tolist()))
                    elif self.data_type.strip() == 'pre':
                        pDlist.append(PREData(atom_type.strip(), resi.strip(),
                        exp, tol, atom.get_coord().tolist()))
        self.parsed = pDlist



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
        vvu_scal = (1./((12*math.pi))*10000)
        self.ax           = float(stdin[7])*vvu_scal
        self.rh           = float(stdin[8])*vvu_scal
        self.Ealpha       = float(stdin[9])
        self.Ebeta        = float(stdin[10])
        self.Egamma       = float(stdin[11])



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
        vvu_scal = (1./((12*math.pi))*10000)
        self.ax           = float(stdin[4])*vvu_scal
        self.rh           = float(stdin[5])*vvu_scal
        self.Ealpha       = float(stdin[6])
        self.Ebeta        = float(stdin[7])
        self.Egamma       = float(stdin[8])

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
        pDList = []
        for i in range(0, len(self.dataset)):
            resi, atom_type, exp, tol =  self.dataset[i].split()
            at1, at2 = atom_type[1], atom_type[0]
            for atom in self.model.get_atoms():
                cur_resi      = str(list(atom.get_parent().get_id())[1])
                cur_atype = atom.get_name()
                if (cur_resi.strip() == resi.strip()) and (at1.strip() ==
                 cur_atype.strip()):
                    c1 = (atom.get_coord().tolist())
                if (cur_resi.strip() == resi.strip()) and (at2.strip() ==
                 cur_atype.strip()):
                    c2 = (atom.get_coord().tolist())
                    pDList.append(RDCData(at1.strip(), resi.strip(), exp, tol,
                     c1, at2.strip(),c2))
        self.parsed = pDList


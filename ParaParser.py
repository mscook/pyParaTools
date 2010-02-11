from Bio.PDB   import *
from numpy     import *
from ParaData  import *
from ParaUtils import *


"""Builds data structures from paramagnetic parameters """
#TODO: Investigate using protected data ie _var
#TODO: Test and check test coverage

class ParaParser:
    def __init__(self, stdin):
        """
        A parser. Takes parameters and builds an object for easy retrival and
        manipulation
        @param stdin: The paramagnetic parameters
        @type stdin:  list of length 4.
        """
        self._data_type     = stdin[1]
        self._pdb_fn        = stdin[2]
        self._para_data_fn  = stdin[3]
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self._pdb_fn)
        self._structure     = structure_in
        self._model         = self._structure[0]
        para_data_in        = open(self._para_data_fn).readlines()
        self._dataset       = para_data_in
        self._parsed        = []


    def getDataType(self):
        """
        Returns the type of paramagnetic data parsed
        """
        return self._data_type

    def getPDBFn(self):
        """
        Returns the PDB file parsed
        """
        return self._pdb_fn

    def getParaDataFn(self):
        """
        Returns the paramagnetic data file parsed
        """
        return self._para_data_fn

    def getStructure(self):
        """
        Returns the Bio.PDB structure object
        """
        return self._structure

    def getDataset(self):
        """
        Returns the dataset (post readlines())
        """
        return self._dataset

    def getNumModels(self):
        """
        Returns the number of models contained in the Bio.PDB structure object
        """
        nmodels = 0
        for model in self._structure:
            nmodels = nmodels+1
        return nmodels

    def getParsed(self):
        """
        Returns a list of post parsed ParaData objects
        """
        return self._parsed


    def setPDBFn(self, pdb_name):
        """
        Set/Change the PDB file to be parsed
        @param pdb_name: The name (and path) to the PDB file
        @type pdb_name:  string
        """
        self._pdb_fn = pdb_name

    def setParaDataFn(self, para_name):
        """
        Set/Change the paramagnetic data file parsed
        @param para_name: The name (and path) to the paramagnetic data file
        @type para_name:  string
        """
        self._para_data_fn = para_name

    def setModel(self, num):
        """
        Change the model (effectively the coordinates) associated with the given
        PDB file and paramagnetic dataset
        @param num: Model number (0->1, 1->2)
        @type num:  int
        """
        self._model = self._structure[num]

    def doParse(self):
        """
        Parse the given the imput parameters, build ParaData objects
        allowing for easy retrival and manipulation
        """
        #OPTIMIZE: This method
        pDlist = []
        for i in range(0, len(self._dataset)):
            res, at, exp, tol =  self._dataset[i].split()
            res, at  = res.strip(),at.strip()
            exp, tol = float(exp), float(tol)
            for atom in self._model.get_atoms():
                c_res = str(list(atom.get_parent().get_id())[1]).strip()
                c_at    = atom.get_name().strip()
                if (c_res == res) and (c_at == at):
                    #NOTE: This is possibly a hack...
                    if self._data_type.strip() == 'pcs':
                        pDlist.append(PCSData(at, int(res), exp, tol, atom.get_coord()))
                    elif self._data_type.strip() == 'pre':
                        pDlist.append(PREData(at, int(res), exp, tol, atom.get_coord()))
                    else:
                        print 'Unsupported type...'
        self._parsed = pDlist


    #TODO: Document the following methods
    def getAllXarray(self):
        xA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getCoordx()
            xA[i] = val
        return xA

    def getAllYarray(self):
        yA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getCoordy()
            yA[i] = val
        return yA

    def getAllZarray(self):
        zA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getCoordz()
            zA[i] = val
        return zA

    def getAllMeasarray(self):
        mA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getVal()
            mA[i] = val
        return mA


    def addErrorGuassianMeas(self, sigma, seed=100):
        import random
        random.seed(seed)
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getVal()
            self._parsed[i].setVal(random.gauss(y,sigma))

    def addErrorFlatMeas(self, sigma, seed=100):
        import random
        random.seed(seed)
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getVal()
            up, low = val+sigma, val-sigma
            self._parsed[i].setVal(random.uniform(up,low))


    def pickleParsed(self, outname):
        import pickle
        pf = file(outname, 'w' )
        pickle.dump(self._parsed, pf )

    def unpickleParsed(self, inname):
        import pickle
        pf = file(inname, 'r' )
        pickle.load(self._parsed, pf )



class PCSParser(ParaParser):
    def __init__(self, stdin):
        """
        A PCS parser. Takes parameters and builds an object for easy retrival
        and manipulation
        @param stdin: PCS parameters
        @type stdin:  list of length 12.
        """
        self._data_type     = stdin[1]
        self._pdb_fn        = stdin[2]
        self._para_data_fn  = stdin[3]
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self._pdb_fn)
        self._structure    = structure_in
        self._model        = self._structure[0]
        para_data_in       = open(self._para_data_fn).readlines()
        self._dataset      = para_data_in
        self._metal_loc    = empty((3))
        self._metal_loc[0] = float(stdin[4])
        self._metal_loc[1] = float(stdin[5])
        self._metal_loc[2] = float(stdin[6])
        self._ax           = ToVVU(float(stdin[7]))
        self._rh           = ToVVU(float(stdin[8]))
        self._Ealpha       = float(stdin[9])
        self._Ebeta        = float(stdin[10])
        self._Egamma       = float(stdin[11])


    def getMetalLoc(self):
        """
        Return the numpy array containing the metal position
        """
        return self._metal_loc

    def getMetalLocx(self):
        """
        Return the metal x coordinate
        """
        return self._metal_loc[0]

    def getMetalLocy(self):
        """
        Return the metal y coordinate
        """
        return self._metal_loc[1]

    def getMetalLocz(self):
        """
        Return the metal z coordinate
        """
        return self._metal_loc[2]

    def getAxial(self):
        """
        Return the X-tensor's axial component in SI units
        """
        return FromVVU(self._ax)

    def getRhombic(self):
        """
        Return the X-tensor's rhombic component in SI units
        """
        return FromVVU(self._rh)

    def getAxialVVU(self):
        """
        Return the X-tensor's axial component in VVU units
        """
        return self._ax

    def getRhombicVVU(self):
        """
        Return the X-tensor's rhombic component in VVU units
        """
        return self._rh

    def getAlpha(self):
        """
        Return the X-tensor frame Alpha angle
        """
        return self._Ealpha

    def getBeta(self):
        """
        Return the X-tensor frame Beta angle
        """
        return self._Ebeta

    def getGamma(self):
        """
        Return the X-tensor frame Gamma angle
        """
        return self._Egamma


    def setMetalLoc(self, metal_xyz):
        """
        Set/Change the metal position
        @param metal_xyz: The metal ion x,y,z coordinates
        @type metal_xyz:  numpy array
        """
        self._metal_loc = metal_xyz

    def setMetalLocx(self, mx):
        """
        Set/Change the metal x coordinate
        @param mx: The x coordinate
        @type mx:  float
        """
        self._metal_loc[0] = mx

    def setMetalLocy(self, my):
        """
        Set/Change the metal y coordinate
        @param my: The y coordinate
        @type my:  float
        """
        self._metal_loc[1] = my

    def setMetalLocz(self, mz):
        """
        Set/Change the metal z coordinate
        @param mz: The x coordinate
        @type mz:  float
        """
        self._metal_loc[2] = mz

    def setAxial(self, axial):
        """
        Set/Change the X-tensor's axial component in VVU units
        @param axial: The axial component in SI units
        @type axial:  float
        """
        self._ax = ToVVU(axial)

    def setRhombic(self, rhombic):
        """
        Set/Change the X-tensor's rhombic component in VVU units
        @param rhombic: The rhombic component in SI units
        @type rhombic:  float
        """
        self._rh = ToVVU(rhombic)

    def setAlpha(self, alpha):
        """
        Set/Change the X-tensor frame Alpha angle
        @param alpha: The Alpha angle in degrees
        @type alpha: float
        """
        self._Ealpha = alpha

    def setBeta(self, beta):
        """
        Set/Change the X-tensor frame Beta angle
        @param beta: The Beta angle in degrees
        @type beta: float
        """
        self._Ebeta = beta

    def setGamma(self, gamma):
        """
        Set/Change the X-tensor frame Gamma angle
        @param gamma: The Gamma angle in degrees
        @type gamma: float
        """
        self._Egamma = gamma



class PREParser(ParaParser):
    def __init__(self, stdin):
        """
        A PRE parser. Takes parameters and builds an object for easy retrival
        and manipulation
        @param stdin: PRE parameters
        @type stdin:  list of length 8.
        """
        self._data_type     = stdin[1]
        self._pdb_fn        = stdin[2]
        self._para_data_fn  = stdin[3]
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self._pdb_fn)
        self._structure    = structure_in
        self._model        = self._structure[0]
        para_data_in       = open(self._para_data_fn).readlines()
        self._dataset      = para_data_in
        self._metal_loc    = zeros((3))
        self._metal_loc[0] = float(stdin[4])
        self._metal_loc[1] = float(stdin[5])
        self._metal_loc[2] = float(stdin[6])
        self._c            = zeros((1))
        self._c[0]         = float(stdin[7])


    def getMetalLoc(self):
        """
        Return the numpy array containing the metal position
        """
        return self._metal_loc

    def getMetalLocx(self):
        """
        Return the metal x coordinate
        """
        return self._metal_loc[0]

    def getMetalLocy(self):
        """
        Return the metal y coordinate
        """
        return self._metal_loc[1]

    def getMetalLocz(self):
        """
        Return the metal z coordinate
        """
        return self._metal_loc[2]

    def getConstant(self):
        """
        Return the value of the PRE constant
        """
        return self._c


    def setMetalLoc(self, metal_xyz):
        """
        Set/Change the metal position
        @param metal_xyz: The metal ion x,y,z coordinates
        @type metal_xyz:  numpy array
        """
        self._metal_loc = metal_xyz

    def setMetalLocx(self, mx):
        """
        Set/Change the metal x coordinate
        @param mx: The x coordinate
        @type mx:  float
        """
        self._metal_loc[0] = mx

    def setMetalLocy(self, my):
        """
        Set/Change the metal y coordinate
        @param my: The y coordinate
        @type my:  float
        """
        self._metal_loc[1] = my

    def setMetalLocz(self, mz):
        """
        Set/Change the metal z coordinate
        @param mz: The z coordinate
        @type mz:  float
        """
        self._metal_loc[2] = mz

    def setConstant(self, constant):
        """
        Set/Change the value of the PRE constant
        @param constant: The PRE constant
        @type constant:  float
        """
        self._c = constant



class RDCParser(ParaParser):
    def __init__(self, stdin):
        """
        A RDC parser. Takes parameters and builds an object for easy retrival
        and manipulation
        @param stdin: RDC parameters
        @type stdin:  list of length 9.
        """
        self._data_type     = stdin[1]
        self._pdb_fn        = stdin[2]
        self._para_data_fn  = stdin[3]
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self._pdb_fn)
        self._structure    = structure_in
        self._model        = self._structure[0]
        para_data_in       = open(self._para_data_fn).readlines()
        self._dataset      = para_data_in
        #TODO: Need to think about what units the alignment tensor is stored in
        self._ax           = ToVVU(float(stdin[4]))
        self._rh           = ToVVU(float(stdin[5]))
        self._Ealpha       = float(stdin[6])
        self._Ebeta        = float(stdin[7])
        self._Egamma       = float(stdin[8])


    def getAxial(self):
        #FIXME: RE: units of the alignment-tensor
        """
        Return the alignment-tensor's axial component in SI units
        """
        return FromVVU(self._ax)

    def getRhombic(self):
        #FIXME: RE: units of the alignment-tensor
        """
        Return the alignment-tensor's rhombic component in SI units
        """
        return FromVVU(self._rh)

    def getAlpha(self):
        """
        Return the alignment-tensor frame Alpha angle
        """
        return self._Ealpha

    def getBeta(self):
        """
        Return the alignment-tensor frame Beta angle
        """
        return self._Ebeta

    def getGamma(self):
        """
        Return the alignment-tensor frame Gamma angle
        """
        return self._Egamma


    def setAxial(self, axial):
        #FIXME: RE: units of the alignment-tensor
        """
        Set/Change the alignment-tensor's axial component in VVU units
        @param axial: The axial component in SI units
        @type axial:  float
        """
        self._ax = ToVVU(axial)

    def setRhombic(self, rhombic):
        #FIXME: RE: units of the alignment-tensor
        """
        Set/Change the alignment-tensor's rhombic component in VVU units
        @param rhombic: The axial component in SI units
        @type rhombic:  float
        """
        self._rh = ToVVU(rhombic)

    def setAlpha(self, alpha):
        """
        Set/Change the alignment-tensor frame Alpha angle
        @param alpha: The Alpha angle in degrees
        @type alpha:  float
        """
        self._Ealpha = alpha

    def setBeta(self, beta):
        """
        Set/Change the alignment-tensor frame Beta angle
        @param beta: The Beta angle in degrees
        @type beta:  float
        """
        self._Ebeta = beta

    def setGamma(self, gamma):
        """
        Set/Change the alignment-tensor frame Gamma angle
        @param gamma: The Alpha angle in degrees
        @type gamma:  float
        """
        self._Egamma = gamma

    def doParse(self):
        """
        Parse the given the imput parameters, build ParaData objects
        allowing for easy retrival and manipulation
        """
        #OPTIMIZE: This method
        pDList = []
        for i in range(0, len(self._dataset)):
            res, at, exp, tol =  self._dataset[i].split()
            at1, at2 = atom_type[1], atom_type[0]
            res = res.strip()
            exp, tol = float(exp), float(tol)
            for atom in self._model.get_atoms():
                c_res      = str(list(atom.get_parent().get_id())[1]).strip()
                c_at = atom.get_name().strip()
                if (c_res == res) and (at1 == c_at):
                    c1 = atom.get_coord()
                if (c_res == res) and (at2 == c_at):
                    c2 = atom.get_coord()
                    pDList.append(RDCData(at1, int(res), exp, tol, c1, at2 ,c2))
        self._parsed = pDList


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
        Change the model (effectively the coordinates) associated with the
        given PDB file and paramagnetic dataset
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
                        pDlist.append(PCSData(at, int(res), exp, tol, \
                                                        atom.get_coord()))
                    elif self._data_type.strip() == 'pre':
                        pDlist.append(PREData(at, int(res), exp, tol, \
                                                        atom.get_coord()))
                    else:
                        print 'Unsupported type...'
        self._parsed = pDlist


    def writeDataSet(self, fname):
        """
        Write the data to file for later use.
        @param fname: The output filename
        """
        #TODO: Check if this method works with RDC.
        fout = open(fname, 'w')
        objsL = self.getParsed()
        for obj in range(0, len(objsL)):
            out = str(objsL[obj].getId()) +'\t'  + str(objsL[obj].getName()) \
            +'\t'+str(objsL[obj].getVal()) +'\t' + str(objsL[obj].getTol())
            fout.write(out+'\n')
        print "Have written data to "+fname


    def getAllXarray(self):
        """
        Returns all the parsed x coordinates in a single array. Exploited in
        fitting.
        """
        xA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getCoordx()
            xA[i] = val
        return xA

    def getAllYarray(self):
        """
        Returns all the parsed y coordinates in a single array. Exploited in
        fitting.
        """
        yA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getCoordy()
            yA[i] = val
        return yA

    def getAllZarray(self):
        """
        Returns all the parsed x coordinates in a single array. Exploited in
        fitting.
        """
        zA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getCoordz()
            zA[i] = val
        return zA

    def getAllMeasarray(self):
        """
        Returns all the parsed measured experimental data in a single array.
        Exploited in fitting.
        """
        mA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getVal()
            mA[i] = val
        return mA

    def getAllTolarray(self):
        """
        Returns all the parsed experimental tolerances in a single array.
        Exploited in fitting.
        """
        tolA = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getTol()
            tolA[i] = val
        return tolA


    def addErrorGuassianMeas(self, delta=1):
        """
        Add random Guassian error to the measured experimental value based on
        the experimental tolerance.
        An optional parameter sigma (delta=1) can be used to adjust sigma
        (defaults to the experimental tolerance)
        """
        import random
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getVal()
            tol = self._parsed[i].getTol()
            print 79*'-'
            print "WARNING: sigma is by default the exprimental tolerance"
            print "When using Guassian ~ 32% of the time the sigma will be"
            print "larger than this tolerance"
            print 79*'-'
            sigma = tol*delta
            self._parsed[i].setVal(random.gauss(val,sigma))

    def addErrorFlatMeas(self, delta=0.01):
        """
        Add error from a "square" distribution to the measured experimental
        value based on relative error (does not account for the experimental
        tolerance.
        An optional parameter sigma (delta=0.01) can be used to adjust the
        relative error.
        """
        import random
        for i in range (0, len(self._parsed)):
            val = self._parsed[i].getVal()
            #print 79*'-'
            #print "WARNING: This method does not consider exprimental tolerance"
            #print "It works best for a 'flat' relative error"
            #print 79*'-'
            sigma = val*delta
            up, low = val+sigma, val-sigma
            self._parsed[i].setVal(random.uniform(up,low))


    def pickleParsed(self, outname):
        """
        Pickle (store the entire parsed datastructure).Mainly useful when
        parsing large PDB files.
        @param outname: the filename to store the pickled datastructure
        @type  outname: string
        """
        import pickle
        pf = file(outname, 'w' )
        pickle.dump(self._parsed, pf )

    def unpickleParsed(self, inname):
        """
        Unpickle (retrive the entire parsed datastructure).Mainly useful when
        parsing large PDB files.
        @param inname: the filename to retrive the pickled object
        @type  inname: string
        """
        import pickle
        pf = file(inname, 'r' )
        pickle.load(self._parsed, pf )

    def setCalcedToObs(self):
        """
        To update a calculated value to become an experimental value
        """
        for i in range (0, len(self._parsed)):
            self._parsed[i].setVal(self._parsed[i].getCVal)



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


    def add2ndSite(self):
        #NOTE: No longer nessacary
        """
        This is a bit of a hack, but in the case of fitting two paramagnetic
        centres we can store info on the second site in the object. The
        following 3 methods (inclusive) deal with this case.
        Initializes a empty datastore
        """
        self._site2 = zeros(8)

    def populate2ndSite(self, xm2,ym2,zm2, ax2,rh2, a2,b2,g2):
        #NOTE: No longer nessacary
        """
        Add X-tensor params to the datastore.
        """
        self._site2[0] = xm2
        self._site2[1] = ym2
        self._site2[2] = zm2
        self._site2[3] = ToVVU(ax2)
        self._site2[4] = ToVVU(rh2)
        self._site2[5] = a2
        self._site2[6] = b2
        self._site2[7] = g2

    def get2ndSite(self):
        #NOTE: No longer nessacary
        """
        Return the X-tensor params from the datastore.
        """
        site2_data = zeros(8)
        site2_data[0] = self._site2[0]
        site2_data[1] = self._site2[1]
        site2_data[2] = self._site2[2]
        site2_data[3] = FromVVU(self._site2[3])
        site2_data[4] = FromVVU(self._site2[4])
        site2_data[5] = self._site2[5]
        site2_data[6] = self._site2[6]
        site2_data[7] = self._site2[7]
        return site2_data


    def getTensorParams(self):
        """
        Return all 8 X-tensor parameters.
        """
        Xt = zeros(8)
        Xt[0] = self._metal_loc[0]
        Xt[1] = self._metal_loc[1]
        Xt[2] = self._metal_loc[2]
        Xt[3] = FromVVU(self._ax)
        Xt[4] = FromVVU(self._rh)
        Xt[5] = self._Ealpha
        Xt[6] = self._Ebeta
        Xt[7] = self._Egamma
        return Xt


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


    def add2ndSite(self):
        #NOTE: No longer nessacary
        """
        This is a bit of a hack, but in the case of fitting two paramagnetic
        centres we can store info on the second site in the object. The
        following 3 methods (inclusive) deal with this case.
        Initializes a empty datastore
        """
        self._site2 = zeros(4)

    def populate2ndSite(self, xm2,ym2,zm2, c2):
        #NOTE: No longer nessacary
        """
        Add X-tensor params to the datastore.
        """
        self._site2[0] = xm2
        self._site2[1] = ym2
        self._site2[2] = zm2
        self._site2[3] = c2

    def get2ndSite(self):
        #NOTE: No longer nessacary
        """
        Return the PRE params from the datastore.
        """
        return self._site2


    def getSiteParams(self):
        """
        Return all 4 PRE centre parameters
        """
        Ps = zeros(4)
        Ps[0] = self._metal_loc[0]
        Ps[1] = self._metal_loc[1]
        Ps[2] = self._metal_loc[2]
        Ps[3] = self._c
        return Ps


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
        self._B0           = 18.792343
        self._temp         = 298.0
        self._order        = 1.0


    def getTensorParams(self):
        """
        Return all 5 alignment-tensor parameters.
        """
        At = zeros(5)
        At[0] = FromVVU(self._ax)
        At[1] = FromVVU(self._rh)
        At[2] = self._Ealpha
        At[3] = self._Ebeta
        At[4] = self._Egamma
        return At

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

    def getB0(self):
        """
        Return the magnetic field strength (in Tesla)
        """
        return self._B0

    def getTemp(self):
        """
        Return the temperature (in Kelvin)
        """
        return self._temp

    def getOrder(self):
        """
        Return the order paramameter
        """
        return self._order

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

    def setB0(self, b0):
        """
        Set the magnetic field strength (in Tesla)
        """
        self._B0 = b0

    def setTemp(self, temperature):
        """
        Set the temperature (in Kelvin)
        """
        self._temp = temperature

    def setOrder(self, Odparam):
        """
        Set the order paramameter
        """
        self._order = Odparam


    def getAllXarray(self):
        """
        Returns all the parsed x (x2) coordinates in a single array (x2).
        Exploited in fitting.
        """
        xA1 = zeros(len(self._parsed))
        xA2 = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val1, val2 = self._parsed[i].getCoordx()
            xA1[i] = val1
            xA2[i] = val2
        return xA1, xA2

    def getAllYarray(self):
        """
        Returns all the parsed y (x2) coordinates in a single array (x2).
        Exploited in fitting.
        """
        yA1 = zeros(len(self._parsed))
        yA2 = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val1, val2 = self._parsed[i].getCoordy()
            yA1[i] = val1
            yA2[i] = val2
        return yA1, yA2

    def getAllZarray(self):
        """
        Returns all the parsed z (x2) coordinates in a single array (x2).
        Exploited in fitting.
        """
        zA1 = zeros(len(self._parsed))
        zA2 = zeros(len(self._parsed))
        for i in range (0, len(self._parsed)):
            val1, val2 = self._parsed[i].getCoordz()
            zA1[i] = val1
            zA2[i] = val2
        return zA1, zA2



    def doParse(self):
        """
        Parse the given the imput parameters, build ParaData objects
        allowing for easy retrival and manipulation
        """
        #OPTIMIZE: This method
        pDList = []
        for i in range(0, len(self._dataset)):
            res, at, exp, tol =  self._dataset[i].split()
            at1, at2 = at[1], at[0]
            #NOTE: This is a hack. Find something more elegant
            at2 = at2+'N'
            res = res.strip()
            exp, tol = float(exp), float(tol)
            for atom in self._model.get_atoms():
                c_res      = str(list(atom.get_parent().get_id())[1]).strip()
                c_at = atom.get_name().strip()
                if (c_res == res) and (at1 == c_at):
                    c1 = atom.get_coord()
                if (c_res == res) and (at2 == c_at):
                    c2 = atom.get_coord()
                    pDList.append(RDCData(at1, int(res), exp, tol, c1, at2[0] ,c2))
        self._parsed = pDList


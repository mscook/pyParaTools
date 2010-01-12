import sys
from Bio.PDB import *

class ParaParser:
    def __init__(self, stdin):
        self.data_type    = stdin[1]
        self.pdb_fn       = stdin[2]
        self.para_data_fn = stdin[3]

    def getDataType(self):
        return self.data_type

    def getPDBFn(self):
        return self.pdb_fn

    def getParaDataFn(self):
        return self.para_data_fn


    def setParaDataFn(self, data_name):
        self.data_type = data_name

    def setPDBFn(self, pdb_name):
        self.pdb_fn = pdb_name

    def setParaDataFn(self, para_name):
        para_data_fn = para_name


    def open_data()
        parser=PDBParser()
        structure_in=parser.get_structure('cur_structure', self.pdb_fn)
        para_data_in = open(self.para_data_fn.readlines())
        return structure_in, para_data_in



    def check_models()
        structure, data = open_data()
        nmodels = 0
        for model in structure:
            nmodels = nmodels+1
        if nmodels > 1:
            print "WARNING:", nmodels,
                "models exist. Working with only the first model"
        print

    def match_data(

        parsed = {}
        for i in range(0, len(data)):
            resi, atom_type, exp, tol =  data[i].split()
            for atom in structure[0].get_atoms():
                cur_resi      = str(list(atom.get_parent().get_id())[1])
                cur_atom_type = atom.get_name()
                if (cur_resi.strip() == resi.strip()) and (cur_atom_type.strip() == atom_type.strip()):
                    d_key =  resi
                    d_val = atom.get_coord().tolist()
                    d_val.insert(0,float(tol))
                    d_val.insert(0,float(exp))
                    d_val.insert(0,atom_type)
                    d_val.insert(0,int(resi))
                results[d_key] = d_val
        #See:http://pyfaq.infogami.com/how-do-i-create-a-multidimensional-list
        return parsed



class PCSParser(ParaParser):
    def __init__(self, stdin):
        self.data_type    = stdin[1]
        self.pdb_fn       = stdin[2]
        self.para_data_fn = stdin[3]
        self.metal_loc    = []
        self.metal_loc.append(float(sys.argv[4])
        self.metal_loc.append(float(sys.argv[5])
        self.metal_loc.append(float(sys.argv[6])
        self.ax           = stdin[7]
        self.rh           = stdin[8]
        self.Ealpha       = stdin[9]
        self.Ebeta        = stdin[10]
        self.Egamma       = stdin[11]

            def setMetalLoc(self, metal_xyz):
        self.metal_loc = []
        self.metal_loc.append(metal_xyz[0])
        self.metal_loc.append(metal_xyz[1])
        self.metal_loc.append(metal_xyz[2])

    def getMetalLoc(self):
        return self.metal_loc

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
    def setAlpha(self, A):
        self.Ealpha = A
    def setBeta(self, B):
        self.Ebeta  = B
    def setGamma(self, G):
        self.Egamma = G


class PREParser(ParaParser):
    def __init__(self, stdin):
        self.data_type    = stdin[1]
        self.pdb_fn       = stdin[2]
        self.para_data_fn = stdin[3]
        self.metal_loc    = []
        self.metal_loc.append(float(sys.argv[4])
        self.metal_loc.append(float(sys.argv[5])
        self.metal_loc.append(float(sys.argv[6])
        self.c            = stdin[7]

        def getConstant(self):
        return self.c

        def setConstant(self, constant):
        self.c = constant



        parser=PDBParser()	# Could use: PERMISSIVE=1 if required
	structure=parser.get_structure('cur_structure', sname)
	data = open(dname, 'r').readlines()

	nmodels = 0
        for model in structure:
                nmodels = nmodels+1
	if nmodels > 1:
		print "WARNING:", nmodels, "models exist. Working with only the first model"
		print

	for i in range(0, len(data)):
		resi, atom_type, exp, tol =  data[i].split()
		for atom in structure[0].get_atoms():
			cur_resi      = str(list(atom.get_parent().get_id())[1])
			cur_atom_type = atom.get_name()
			if (cur_resi.strip() == resi.strip()) and (cur_atom_type.strip() == atom_type.strip()):
				d_key =  resi
				d_val = atom.get_coord().tolist()
				d_val.insert(0,float(tol))
				d_val.insert(0,float(exp))
				d_val.insert(0,atom_type)
				d_val.insert(0,int(resi))
				results[d_key] = d_val
	return results


class PCSParser(ParaParser)

class PRE



type = [1]
pdb     = sys.argv[2]
data    = sys.argv[3]
metal   = [float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])]
ax      = float(sys.argv[7])
rh      = float(sys.argv[8])

euler   = [float(sys.argv[9]), float(sys.argv[10]), float(sys.argv[11])]
euler[0]= math.radians(euler[0])
euler[1]= math.radians(euler[1])
euler[2]= math.radians(euler[2])


    def getName(self):
        return self.resi_name
    def getId(self):
        return self.resi_id
    def getVal(self):
        return self.exp_val
    def getTol(self):
        return self.e_tol
    def getCoord(self):
        return self.s_coord
    def getType(self):
        return self.p_type

    def setName(self, name):
        self.resi_name = name
    def setId(self, r_id):
        self.resi_id = r_id
    def setVal(self, val):
        self.exp_val = val
    def setTol(self, tol):
        self.e_tol = tol
    def setCoord(self, coord):
        self.s_coord = coord


class ParaData:

    def __init__(self, resi_name, resi_id, exp_val, e_tol, s_coord):
        self.resi_name = resi_name
        self.resi_id   = resi_id
        self.exp_val   = exp_val
        self.e_tol     = e_tol
        self.s_coord   = s_coord
        self.p_type    = 'unknown'
        
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


class PCSData(ParaData):
    def __init__(self, resi_name, resi_id, exp_val, e_tol, s_coord):
        self.resi_name  = resi_name
        self.resi_id    = resi_id
        self.exp_val    = exp_val
        self.e_tol      = e_tol
        self.s_coord    = s_coord
        self.p_type     = 'pcs'
        
        
class PREData(ParaData):
    def __init__(self, resi_name, resi_id, exp_val, e_tol, s_coord):
        self.resi_name  = resi_name
        self.resi_id    = resi_id
        self.exp_val    = exp_val
        self.e_tol      = e_tol
        self.s_coord    = s_coord
        self.p_type     = 'pre'
        
        
class RDCData(ParaData):
    def __init__(self, resi_name, resi_name2, resi_id, exp_val, e_tol, s_coord, s_coord2):
        self.resi_name  = resi_name
        self.resi_name2 = resi_name2
        self.resi_id    = resi_id
        self.exp_val    = exp_val
        self.e_tol      = e_tol
        self.s_coord    = s_coord
        self.s_coord2   = s_coord2
        self.p_type     = 'rdc'
        
    def getName(self):
        name = self.resi_name+"_"+self.resi_name2
        return name
    def getCoord(self):
        return self.s_coord, self.s_coord2

    def setName2(self, name2):
        self.resi_name2 = name2
    def setCoord2(self, coord2):
        self.s_coord2 = coord2
        
    def lookupMGR(self):
        #Need to add 
        #mgr = {'C':'', 'N':19.338e6, 'H':267.522e6}
        
class ParaStore:
...     def __init__(self, initial_student_list):
...         self.student_names = initial_student_list
...         self.students = {}
...         for name in self.student_names:
...             self.students[name] = Student(name)
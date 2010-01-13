class ParaData:

    def __init__(self, spin_type, spin_id, exp_val, e_tol, s_coord):
        """
        A paramagnetic data type
        @param spin_type: The spin type
        @type spin_type:  string
        @param spin_id:   The residue id
        @type spin_id:    int
        @param exp_val:   The experimental value
        @type exp_val:    float
        @param e_tol:     The experimental tolerance
        @type e_tol:      float
        @param s_coord:   The coordinates of the spin
        @type s_coord:    list
        """
        self.spin_type = spin_type
        self.spin_id   = spin_id
        self.exp_val   = exp_val
        self.e_tol     = e_tol
        self.s_coord   = s_coord
        self.p_type    = 'unknown'

    def getName(self):
        """
        Returns the spin type assigned to this spin
        """
        return self.spin_type

    def getId(self):
        """
        Returns the residue number assigned to this spin
        """
        return self.spin_id

    def getVal(self):
        """
        Returns the experimental value assigned to this spin
        """
        return self.exp_val

    def getTol(self):
        """
        Returns the experimental tolerance assigned to this spin
        """
        return self.e_tol

    def getCoord(self):
        """
        Returns the coordinates (x,y,z) assigned to this spin
        """
        return self.s_coord

    def getCoordx(self):
        """
        Returns the x coordinate assigned to this spin
        """
        return self.s_coord[0]

    def getCoordy(self):
        """
        Returns the y coordinate assigned to this spin
        """
        return self.s_coord[1]

    def getCoordz(self):
        """
        Returns the z coordinate assigned to this spin
        """
        return self.s_coord[2]

    def getType(self):
        """
        Returns the experimental type (PCS, PRE, RDC) assigned to this spin
        """
        return self.p_type


    def setName(self, name):
        """
        Update the spin type
        @param name: The atom type
        @type name:  string
        """
        self.spin_type = name

    def setId(self, r_id):
        """
        Update the residue id
        @param r_id: The residue id
        @type r_id: int
        """
        self.spin_id = r_id

    def setVal(self, val):
        """
        Update the experimental value
        @param val: The experimental value
        @type val: float
        """
        self.exp_val = val

    def setTol(self, tol):
        """
        Update the experimental tolerance
        @param tol: The experimental tolerance
        @type tol: float
        """
        self.e_tol = tol

    def setCoord(self, coord):
        """
        Update the coordinates of the spin
        @param coord: The coordinates of the spin
        @type coord: list
        """
        self.s_coord = coord

    def setCoordx(self, coord_x):
        """
        Update the x coordinate of the spin
        @param coord: The x coordinate of the spin
        @type coord: float
        """
        self.s_coord[0] = coord_x

    def setCoordy(self, coord_y):
        """
        Update the y coordinate of the spin
        @param coord: The y coordinate of the spin
        @type coord: float
        """
        self.s_coord[1] = coord_y

    def setCoordz(self, coord_z):
        """
        Update the z coordinate of the spin
        @param coord: The z coordinate of the spin
        @type coord: float
        """
        self.s_coord[2] = coord_z

    def printNice(self):
        """
        Print 'Numbat' like format
        """
        print '%s%6s%8.3f%8.3f' % (self.spin_type, self.spin_id,
        self.exp_val, self.e_tol)


class PCSData(ParaData):
    def __init__(self, spin_type, spin_id, exp_val, e_tol, s_coord):
        """
        A PCS data container
        @param spin_type: The spin type
        @type spin_type:  string
        @param spin_id:   The residue id
        @type spin_id:    int
        @param exp_val:   The experimental PCS value
        @type exp_val:    float
        @param e_tol:     The experimental PCS tolerance
        @type e_tol:      float
        @param s_coord:   The coordinates of the spin
        @type s_coord:    list
        """
        self.spin_type  = spin_type
        self.spin_id    = spin_id
        self.exp_val    = exp_val
        self.e_tol      = e_tol
        self.s_coord    = s_coord
        self.p_type     = 'pcs'



class PREData(ParaData):
    def __init__(self, spin_type, spin_id, exp_val, e_tol, s_coord):
        """
        A PRE data container
        @param spin_type: The spin type
        @type spin_type:  string
        @param spin_id:   The residue id
        @type spin_id:    int
        @param exp_val:   The experimental PRE value
        @type exp_val:    float
        @param e_tol:     The experimental PRE tolerance
        @type e_tol:      float
        @param s_coord:   The coordinates of the spin
        @type s_coord:    list
        """
        self.spin_type  = spin_type
        self.spin_id    = spin_id
        self.exp_val    = exp_val
        self.e_tol      = e_tol
        self.s_coord    = s_coord
        self.p_type     = 'pre'



class RDCData(ParaData):
    def __init__(self, spin_type, spin_id, exp_val, e_tol,
     s_coord, spin_type2, s_coord2):
        """
        A RDC data container
        @param spin_type: The spin type of the 1st atom in the coupling
        @type spin_type:  string
        @param spin_type2:The spin type of the 2nd atom in the coupling
        @type spin_type2: string
        @param spin_id:   The residue id
        @type spin_id:    int
        @param exp_val:   The experimental PCS value
        @type exp_val:    float
        @param e_tol:     The experimental PCS tolerance
        @type e_tol:      float
        @param s_coord:   The coordinates of the first spin in the coupling
        @type s_coord:    list
        @param s_coord2:  The coordinates of the second spin in the coupling
        @type s_coord2:   list
        """
        self.spin_type  = spin_type
        self.spin_id    = spin_id
        self.exp_val    = exp_val
        self.e_tol      = e_tol
        self.s_coord    = s_coord
        self.spin_type2 = spin_type2
        self.s_coord2   = s_coord2
        self.p_type     = 'rdc'

    def getName(self):
        """
        Return the coupling type
        """
        name = self.spin_type+self.spin_type2
        return name

    def getCoord(self):
        """
        Returns the coordinates of *both* spins in the coupling
        """
        return self.s_coord, self.s_coord2

    def getCoordx(self):
        """
        Returns the x coordinates of *both* spins in the coupling
        """
        return self.s_coord[0], self.s_coord2[0]

    def getCoordy(self):
        """
        Returns the y coordinates of *both* spins in the coupling
        """
        return self.s_coord[1], self.s_coord2[1]

    def getCoordz(self):
        """
        Returns the z coordinates of *both* spins in the coupling
        """
        return self.s_coord[2], self.s_coord2[2]


    def setName2(self, name2):
        """
        Update the spin type for the second spin
        @param name2: The atom type
        @type name2: string
        """
        self.spin_type2 = name2

    def setCoord2(self, coord2):
        """
        Update the coordinates for the second spin
        @param coord2: The spin cordinates
        @type coord2: list
        """
        self.s_coord2 = coord2

    def setCoord2x(self, coord2_x):
        """
        Update the x coordinate for the second spin
        @param coord2: The cordinate
        @type coord2: float
        """
        self.s_coord2 = coord2_x

    def setCoord2y(self, coord2_y):
        """
        Update the y coordinate for the second spin
        @param coord2: The cordinate
        @type coord2: float
        """
        self.s_coord2 = coord2_y

    def setCoord2z(self, coord2_z):
        """
        Update the z coordinate for the second spin
        @param coord2: The cordinate
        @type coord2: float
        """
        self.s_coord2 = coord2_z

    def printNice(self):
        """
        Print 'Numbat' like format
        """
        print '%s%6s%8.3f%8.3f' % (self.spin_type+self.spin_type2, self.spin_id,
        self.exp_val, self.e_tol)


    def lookupMGR(self):
        """
        Return the gyromagnetic ratios for the coupling.
        See: http://nmrwiki.org/wiki/index.php?title=Gyromagnetic_ratio
        """
        PI2        = 2*3.1415926535897931
        H1mgr     = (PI2*42.576)*1e6
        C13mgr    = (PI2*10.705)*1e6
        Nmgr = []
        N14mgr    = (PI2*3.0766)*1e6
        N15mgr    = (PI2*-4.315)*1e6
        Nmgr.append(N14mgr)
        Nmgr.append(N15mgr)
        O17mgr    = (PI2*-5.7716)*1e6
        mgr = {'H':H1mgr, 'C':C13mgr, 'N':Nmgr, 'O':O17mgr}
        mgr_1 = mgr[self.spin_type]
        mgr_2 = mgr[self.spin_type2]
        return mgr_1, mgr_2


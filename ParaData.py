from     numpy import *
#TODO: Check above is nesacarry

"""A paramagnetic datatype friendly container """
    #TODO: Investigate using protected data ie _var - *DONE*
    #TODO: Test and check current test coverage

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
        @type s_coord:    numpy array
        """
        self._spin_type = spin_type
        self._spin_id   = spin_id
        self._exp_val   = exp_val
        self._calced    = 0.0
        self._e_tol     = e_tol
        self._s_coord   = s_coord
        self._p_type    = 'unknown'


#    def __repr__(self):
#        #TODO: Check that this method actually functions
#        return '%s%12i%12.3f%12.3f\n' % (self._spin_type, self._spin_id,
#        self._exp_val, self._e_tol)


    def getName(self):
        """
        Returns the spin type assigned to this spin
        """
        return self._spin_type

    def getId(self):
        """
        Returns the residue number assigned to this spin
        """
        return self._spin_id

    def getVal(self):
        """
        Returns the experimental value assigned to this spin
        """
        return self._exp_val

    def getCVal(self):
        """
        Returns the calculated value assigned to this spin
        """
        return self._calced

    def getTol(self):
        """
        Returns the experimental tolerance assigned to this spin
        """
        return self._e_tol

    def getCoord(self):
        """
        Returns the coordinates (x,y,z) assigned to this spin
        """
        return self._s_coord

    def getCoordx(self):
        """
        Returns the x coordinate assigned to this spin
        """
        return self._s_coord[0]

    def getCoordy(self):
        """
        Returns the y coordinate assigned to this spin
        """
        return self._s_coord[1]

    def getCoordz(self):
        """
        Returns the z coordinate assigned to this spin
        """
        return self._s_coord[2]

    def getType(self):
        """
        Returns the experimental type (PCS, PRE, RDC) assigned to this spin
        """
        return self._p_type


    def setName(self, name):
        """
        Update the spin type
        @param name: The atom type
        @type name:  string
        """
        self._spin_type = name

    def setId(self, r_id):
        """
        Update the residue id
        @param r_id: The residue id
        @type r_id:  int
        """
        self._spin_id = r_id

    def setVal(self, val):
        """
        Update the experimental value
        @param val: The experimental value
        @type val:  float
        """
        self._exp_val = val

    def setCVal(self, cval):
        """
        Sets the calculated value assigned to this spin
        """
        self._calced = cval


    def setTol(self, tol):
        """
        Update the experimental tolerance
        @param tol: The experimental tolerance
        @type tol:  float
        """
        self._e_tol = tol

    def setCoord(self, coord):
        """
        Update the coordinates of the spin
        @param coord: The coordinates of the spin
        @type coord:  numpy array
        """
        self._s_coord = coord

    def setCoordx(self, coord_x):
        """
        Update the x coordinate of the spin
        @param coord_x: The x coordinate of the spin
        @type coord_x : float
        """
        self._s_coord[0] = coord_x

    def setCoordy(self, coord_y):
        """
        Update the y coordinate of the spin
        @param coord_y: The y coordinate of the spin
        @type coord_y:  float
        """
        self._s_coord[1] = coord_y

    def setCoordz(self, coord_z):
        """
        Update the z coordinate of the spin
        @param coord_z: The z coordinate of the spin
        @type coord_z:  float
        """
        self._s_coord[2] = coord_z


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
        @type s_coord:    numpy array
        """
        self._spin_type  = spin_type
        self._spin_id    = spin_id
        self._exp_val    = exp_val
        self._calced     = 0.0
        self._e_tol      = e_tol
        self._s_coord    = s_coord
        self._p_type     = 'pcs'



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
        @param calced:    The calculated PRE value
        @type calced:     float
        @param e_tol:     The experimental PRE tolerance
        @type e_tol:      float
        @param s_coord:   The coordinates of the spin
        @type s_coord:    numpy array
        """
        self._spin_type  = spin_type
        self._spin_id    = spin_id
        self._exp_val    = exp_val
        self._calced     = 0.0
        self._e_tol      = e_tol
        self._s_coord    = s_coord
        self._p_type     = 'pre'



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
        @type s_coord:    numpy array
        @param s_coord2:  The coordinates of the second spin in the coupling
        @type s_coord2:   numpy array
        """
        self._spin_type  = spin_type
        self._spin_id    = spin_id
        self._exp_val    = exp_val
        self._calced     = 0.0
        self._e_tol      = e_tol
        self._s_coord    = s_coord
        self._spin_type2 = spin_type2
        self._s_coord2   = s_coord2
        self._p_type     = 'rdc'

    def __repr__(self):
        #TODO: Check that this method actually functions
        return '%s%12i%12.3f%12.3f\n' % (self._spin_type+self._spin_type2,
        self._spin_id, self._exp_val, self._e_tol)


    def getName(self):
        """
        Return the coupling type
        """
        name = self._spin_type+self._spin_type2
        return name

    def getCoord(self):
        """
        Returns the coordinates of *both* spins in the coupling
        """
        return self._s_coord, self._s_coord2

    def getCoordx(self):
        """
        Returns the x coordinates of *both* spins in the coupling
        """
        return self._s_coord[0], self._s_coord2[0]

    def getCoordy(self):
        """
        Returns the y coordinates of *both* spins in the coupling
        """
        return self._s_coord[1], self._s_coord2[1]

    def getCoordz(self):
        """
        Returns the z coordinates of *both* spins in the coupling
        """
        return self._s_coord[2], self._s_coord2[2]


    def setName2(self, name2):
        """
        Update the spin type for the second spin
        @param name2: The atom type
        @type name2: string
        """
        self._spin_type2 = name2

    def setCoord2(self, coord2):
        """
        Update the coordinates for the second spin
        @param coord2: The spin cordinates
        @type coord2:  numpy array
        """
        self._s_coord2 = coord2

    def setCoord2x(self, coord2_x):
        """
        Update the x coordinate for the second spin
        @param coord2_x: The cordinate
        @type coord2_x:  float
        """
        self._s_coord2 = coord2_x

    def setCoord2y(self, coord2_y):
        """
        Update the y coordinate for the second spin
        @param coord2_y: The cordinate
        @type coord2_y:  float
        """
        self._s_coord2 = coord2_y

    def setCoord2z(self, coord2_z):
        """
        Update the z coordinate for the second spin
        @param coord2_z: The cordinate
        @type coord2_z:  float
        """
        self._s_coord2 = coord2_z


"""Utility methods for paramagnetic observables """

import math


def ZXZRot(A, B, G):
    #FIXME: Rewrite to utilize Numpy arrays
    """
     Builds the ZXZ rotation matrix given 3 Euler Angles. See:
         http://mathworld.wolfram.com/EulerAngles.html
     @param A: The (A)lpha angle
     @type A : float
     @param B: The (B)eta angle
     @type B : float
     @param G: The (G)amma angle
     @type G : float
    """
    rot = [None]*3
    for i in range(3):
        rot[i] = [None] * 3
    ca = math.cos(math.radians(A))
    cb = math.cos(math.radians(B))
    cg = math.cos(math.radians(G))
    sa = math.sin(math.radians(A))
    sb = math.sin(math.radians(B))
    sg = math.sin(math.radians(G))
    rot[0][0] = ( cg * ca) - (cb * sa * sg)
    rot[0][1] = ( cg * sa) + (cb * ca * sg)
    rot[0][2] = ( sg * sb)
    rot[1][0] = (-sg * ca) - (cb * sa * cg)
    rot[1][1] = (-sg * sa) + (cb * ca * cg)
    rot[1][2] = ( cg * sb)
    rot[2][0] = ( sb * sa)
    rot[2][1] = (-sb * ca)
    rot[2][2] =   cb
    return rot


def ZYZRot(A, B, G):
    #FIXME: Rewrite to utilize Numpy arrays
    """
    .Builds the ZYZ rotation matrix given 3 Euler Angles. See:
         http://mathworld.wolfram.com/EulerAngles.html
     @param A: The (A)lpha angle
     @type A : float
     @param B: The (B)eta angle
     @type B : float
     @param G: The (G)amma angle
     @type G : float
    """
    rot = [None]*3
    for i in range(3):
        rot[i] = [None] * 3
    ca = math.cos(math.radians(A))
    cb = math.cos(math.radians(B))
    cg = math.cos(math.radians(G))
    sa = math.sin(math.radians(A))
    sb = math.sin(math.radians(B))
    sg = math.sin(math.radians(G))
    rot[0][0] = (-sg * sa) + (cb * ca * cg)
    rot[0][1] = ( sg * ca) + (cb * sa * cg)
    rot[0][2] = ( -cg * sb)
    rot[1][0] = (-cg * sa) - (cb * ca * sg)
    rot[1][1] = (cg * ca) - (cb * sa * sg)
    rot[1][2] = ( sg * sb)
    rot[2][0] = ( sb * ca)
    rot[2][1] = (sb * sa)
    rot[2][2] =   cb
    return rot



def FromVVU(AxorRh):
    """
     Convert from van Vleck Units (vvu = m3/3.77 10-35)
     @param AxorRh: Axial or Rhombic component
     @type AxorRh : float
    """
    return AxorRh/(1./((12*math.pi))*10000)


def ToVVU(AxorRh):
    """
     Convert to van Vleck Units (vvu = m3/3.77 10-35)
     @param AxorRh: Axial or Rhombic component
     @type AxorRh : float
    """
    return AxorRh*(1./((12*math.pi))*10000)


def FixAngle(angle):
    """
     To fix up the angles after optimization as they are not [0:2pi] bound
     @param angle: An Euler angle determined from the optimization
     @type angle:  float
    """
    fix_range = 0.0
    if angle < 0:
       fix_range = 360.0
    return (angle%360)


def lookupMGR(spin_type):
    """
    Return the gyromagnetic ratios for the coupling.
    See: http://nmrwiki.org/wiki/index.php?title=Gyromagnetic_ratio
    """
    #TODO: These need to be checked
    PI2        = 2*math.pi
    H1mgr     = (PI2*42.576)*1e6
    C13mgr    = (PI2*10.705)*1e6
    Nmgr = []
    N14mgr    = (PI2*3.0766)*1e6
    N15mgr    = (PI2*-4.315)*1e6
    Nmgr.append(N14mgr)
    Nmgr.append(N15mgr)
    O17mgr    = (PI2*-5.7716)*1e6
    mgr = {'H':H1mgr, 'C':C13mgr, 'N':Nmgr, 'O':O17mgr}
    return mgr[spin_type]


def getConstants(c_type):
    """
    Lookup NMR constants. For RDC calculations
         @param A: The (A)lpha angle
     @type A : float
     @param B: The (B)eta angle

    """
    #TODO: These need to be checked
    if   c_type == 'hbar':
        return 1.05457148e-34
    elif c_type == 'kboltz':
        return 1.3806503e-23
    elif c_type == distNH:
        return 1.03e-10
    elif c_type == vvuSRDC:
        return 3.76991118431e-35
    else:
        return 0.0


def RDCScal(B0, temp, gH, gN):
    """
    Calculate the RDC scaling factor
    """
    #TODO: I cant remember where/how this works. Check this
    a = -1*(B0**2)*gH*gN*(getConstants(hbar))
    b = ((getConstants(distHN)**3)*getConstants(kboltz)*temp)*(120*math.pi**2)
    return (a/b)*getConstants(vvuSRDC)


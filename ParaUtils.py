"""Utility methods for paramagnetic observables """

import math

def ZXZRot(A, B, G):
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


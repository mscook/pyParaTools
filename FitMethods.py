from     numpy import *
from ParaUtils import *

"""Resual error functions evaluated during the minimization"""

#NOTE: Function name read as: PRE(1M)ODEL(1S)ITE(F)IXED(C)
def PRE1M1SFC(p0, meas, c, x,y,z):
    """
     Optimize for a single PRE centre to a single model with the constant in the
     PRE equation known.
     @param p0:   A numpy array with estimates for the coordinates of the PRE
         centre
     @param meas: A numpy array of measured PRE values
     @param c:    The known/predetermined constant in the PRE equation
     @param x:    A numpy array of x coordinates for the given PRE values
     @param y:    A numpy array of y coordinates for the given PRE values
     @param z:    A numpy array of z coordinates for the given PRE values
    """
    xm,ym,zm = p0
    r_v1 = sqrt((x-xm)**2 +(y-ym)**2 + (z-zm)**2)
    err = meas - c/r_v1**6
    return err

#NOTE: Function name read as: PRE(1M)ODEL(1S)ITE(O)PTIMIZE(C)
def PRE1M1SOC(p0, meas, x,y,z):
    """
     Optimize for a single PRE centre and the constant in the PRE equation to
     a single model. Other parameters are as PRE1M1SFC
     @param p0: A numpy array with estimates for the coordinates of the PRE
         centre and an estimate for the constant in the PRE equations
    """
    xm,ym,zm,c = p0
    r_v1 = sqrt((x-xm)**2 +(y-ym)**2 + (z-zm)**2)
    err = meas - c/r_v1**6
    return err

#NOTE: Function name read as: PRE(2M)ODELS(1S)ITE(F)IXED(C)
def PRE2M1SFC(p0, meas, c, x,y,z, x2,y2,z2):
    """
     Optimize for a single PRE centre to 2 models (dimer) with the constant in
     the PRE equation known. Other parameters are as PRE1M1SFC
     @param x2: A numpy array of 2nd x coordinates for the given PRE values
     @param y2: A numpy array of 2nd y coordinates for the given PRE values
     @param z2: A numpy array of 2nd z coordinates for the given PRE values
    """
    xm,ym,zm = p0
    r_v1 = sqrt( (x-xm)**2 +( y-ym)**2 + ( z-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = meas - (c/r_v1**6 + c/r_v2**6)
    return err

#NOTE: Function name read as: PRE(2M)ODELS(1S)ITE(O)PTIMIZE(C)
def PRE2M1SOC(p0, meas, x,y,z, x2,y2,z2):
    """
     Optimize for a single PRE centre and the constant in the PRE equation to
     2 models. Parameters as PRE1M1SOC and PRE2M1SFC
    """
    xm,ym,zm, c = p0
    r_v1 = sqrt(( x-xm)**2 +( y-ym)**2 + ( z-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = meas - (c/r_v1**6 + c/r_v2**6)
    return err

#NOTE: Function name read as: PRE(1M)ODELS(2S)SITE(F)IXED(C)
def PRE1M2SFC(p0, meas, c, x,y,z):
    """
     Optimize for two PRE centres (dimer case) to a single model with the
     constant in the PRE equation known.
     @param p0:   A numpy array with estimates for the coordinates of each of
        the two PRE centres in the form [<x,y,z>_1  <x,y,z>_2]
    """
    xm, ym, zm, xm2,ym2,zm2 = p0
    r_v1 = sqrt((x -xm )**2 +(y-ym )**2 + (z-zm )**2)
    r_v2 = sqrt((x -xm2)**2 +(y-ym2)**2 + (z-zm2)**2)
    err = meas - (c/r_v1**6 + c/r_v2**6)
    return err

#NOTE: Function name read as: PRE(1M)ODELS(2S)SITE(O)PTIMIZE(C)
def PRE1M2SOC(p0, meas, x,y,z):
    """
     Optimize for two PRE centres (dimer case) to a single model with the
     constant in the PRE equation unknown.
     @param p0:   A numpy array with estimates for the coordinates of each of
        the two PRE centres in the form [<x,y,z>_1  <x,y,z>_2]
    """
    xm, ym, zm, xm2,ym2,zm2, c = p0
    r_v1 = sqrt((x- xm)**2 +(y- ym)**2 + (z- zm)**2)
    r_v2 = sqrt((x-xm2)**2 +(y-ym2)**2 + (z-zm2)**2)
    err = meas - (c/r_v1**6 + c/r_v2**6)
    return err



#NOTE: Function name read as: PCS(1M)ODELS(1S)SITE
def PCS1M1S(p0, meas, x,y,z):
    """
     Optimize for the X-tensor given a single model.
     @param p0: A list of initial estimates for the 8 unknown X-tensor params
    """
    xm, ym, zm, ax, rh, a, b, g = p0
    rot  = ZYZRot(a, b, g)
    X   = x - xm
    Y   = y - ym
    Z   = z - zm
    x_t = rot[0][0]*X + rot[0][1]*Y + rot[0][2]*Z
    y_t = rot[1][0]*X + rot[1][1]*Y + rot[1][2]*Z
    z_t = rot[2][0]*X + rot[2][1]*Y + rot[2][2]*Z
    r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
    r5 = (r2*r2) * sqrt(r2)
    tmp = 1.0/r5
    err = meas - (tmp*(ax * (3.0*z_t*z_t -r2) + rh*1.5*(x_t*x_t - y_t*y_t)))
    return err

#NOTE: Function name read as: PCS(1M)ODELS(1S)SITE(F)IXED(M)ETAL
def PCS1M1SFM(p0, meas, xm,ym,zm, x,y,z):
    """
     Optimize for the X-tensor given a single model with fixed metal position
     @param p0: A list of initial estimates for the 5 unknown X-tensor params
    """
    ax, rh, a, b, g = p0
    rot  = ZYZRot(a, b, g)
    X   = x - xm
    Y   = y - ym
    Z   = z - zm
    x_t = rot[0][0]*X + rot[0][1]*Y + rot[0][2]*Z
    y_t = rot[1][0]*X + rot[1][1]*Y + rot[1][2]*Z
    z_t = rot[2][0]*X + rot[2][1]*Y + rot[2][2]*Z
    r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
    r5 = (r2*r2) * sqrt(r2)
    tmp = 1.0/r5
    err = meas - (tmp*(ax * (3.0*z_t*z_t -r2) + rh*1.5*(x_t*x_t - y_t*y_t)))
    return err

#NOTE: Function name read as: PCS(2M)ODELS(1S)SITE
def PCS2M1S(p0, meas, x,y,z, x2,y2,z2):
    """
     Optimize for the X-tensor given 2 models.
     @param p0: A list of initial estimates for the 8 unknown X-tensor params
    """
    xm, ym, zm, ax, rh, a, b, g = p0
    rot  = ZYZRot(a, b, g)
    X1   = x -  xm
    Y1   = y -  ym
    Z1   = z -  zm
    X2   = x2 - xm
    Y2   = y2 - ym
    Z2   = z2 - zm
    x_t1 = rot[0][0]*X1 + rot[0][1]*Y1 + rot[0][2]*Z1
    y_t1 = rot[1][0]*X1 + rot[1][1]*Y1 + rot[1][2]*Z1
    z_t1 = rot[2][0]*X1 + rot[2][1]*Y1 + rot[2][2]*Z1
    x_t2 = rot[0][0]*X2 + rot[0][1]*Y2 + rot[0][2]*Z2
    y_t2 = rot[1][0]*X2 + rot[1][1]*Y2 + rot[1][2]*Z2
    z_t2 = rot[2][0]*X2 + rot[2][1]*Y2 + rot[2][2]*Z2
    r2_1 = (x_t1*x_t1)+(y_t1*y_t1)+(z_t1*z_t1)
    r2_2 = (x_t2*x_t2)+(y_t2*y_t2)+(z_t2*z_t2)
    r5_1 = (r2_1*r2_1) * sqrt(r2_1)
    r5_2 = (r2_2*r2_2) * sqrt(r2_2)
    tmp_1 = 1.0/r5_1
    tmp_2 = 1.0/r5_2
    pcs1 = (tmp_1*(ax * (3.0*z_t1*z_t1 -r2_1) + rh*1.5*(x_t1*x_t1 - y_t1*y_t1)))
    pcs2 = (tmp_2*(ax * (3.0*z_t2*z_t2 -r2_2) + rh*1.5*(x_t2*x_t2 - y_t2*y_t2)))
    err = meas - (pcs1 + pcs2)
    return err

#NOTE: Function name read as: PCS(2M)ODELS(1S)SITE(F)IXED(M)ETAL
def PCS2M1SFM(p0, meas, xm,ym,zm, x,y,z, x2,y2,z2):
    """
     Optimize for the X-tensor given 2 models with fixed metal position
     @param p0: A list of initial estimates for the 5 unknown X-tensor params
    """
    ax, rh, a, b, g = p0
    rot  = ZYZRot(a, b, g)
    X1   = x -  xm
    Y1   = y -  ym
    Z1   = z -  zm
    X2   = x2 - xm
    Y2   = y2 - ym
    Z2   = z2 - zm
    x_t1 = rot[0][0]*X1 + rot[0][1]*Y1 + rot[0][2]*Z1
    y_t1 = rot[1][0]*X1 + rot[1][1]*Y1 + rot[1][2]*Z1
    z_t1 = rot[2][0]*X1 + rot[2][1]*Y1 + rot[2][2]*Z1
    x_t2 = rot[0][0]*X2 + rot[0][1]*Y2 + rot[0][2]*Z2
    y_t2 = rot[1][0]*X2 + rot[1][1]*Y2 + rot[1][2]*Z2
    z_t2 = rot[2][0]*X2 + rot[2][1]*Y2 + rot[2][2]*Z2
    r2_1 = (x_t1*x_t1)+(y_t1*y_t1)+(z_t1*z_t1)
    r2_2 = (x_t2*x_t2)+(y_t2*y_t2)+(z_t2*z_t2)
    r5_1 = (r2_1*r2_1) * sqrt(r2_1)
    r5_2 = (r2_2*r2_2) * sqrt(r2_2)
    tmp_1 = 1.0/r5_1
    tmp_2 = 1.0/r5_2
    pcs1 = (tmp_1*(ax * (3.0*z_t1*z_t1 -r2_1) + rh*1.5*(x_t1*x_t1 - y_t1*y_t1)))
    pcs2 = (tmp_2*(ax * (3.0*z_t2*z_t2 -r2_2) + rh*1.5*(x_t2*x_t2 - y_t2*y_t2)))
    err = meas - (pcs1 + pcs2)
    return err




def PCSMon2SErr(p, y_v, x1,y1, z1):
    """
     FIXME
     @param p:
     @type p:
     @param y_v:
     @type y_v:
     @param x1:
     @type x1:
     @param y1:
     @type y1:
     @param z1:
     @type z1:
    """
    xm1, ym1, zm1, xm2, ym2, zm2, ax1, rh1, ax2, rh2, a1, b1, g1, a2, b2, g2 = p
    zyz_rot1  = ZYZRot(a1, b1, g1)
    zyz_rot2  = ZYZRot(a2, b2, g2)
    X1, X2   = (x1 - xm1), (x1 - xm2)
    Y1, Y2   = (y1 - ym1), (y1 - ym2)
    Z1, Z2   = (z1 - zm1), (z1 - zm2)
    x_t1 = zyz_rot1[0][0]*X1 + zyz_rot1[0][1]*Y1 + zyz_rot1[0][2]*Z1
    y_t1 = zyz_rot1[1][0]*X1 + zyz_rot1[1][1]*Y1 + zyz_rot1[1][2]*Z1
    z_t1 = zyz_rot1[2][0]*X1 + zyz_rot1[2][1]*Y1 + zyz_rot1[2][2]*Z1
    x_t2 = zyz_rot2[0][0]*X2 + zyz_rot2[0][1]*Y2 + zyz_rot2[0][2]*Z2
    y_t2 = zyz_rot2[1][0]*X2 + zyz_rot2[1][1]*Y2 + zyz_rot2[1][2]*Z2
    z_t2 = zyz_rot2[2][0]*X2 + zyz_rot2[2][1]*Y2 + zyz_rot2[2][2]*Z2
    r2_1 = (x_t1*x_t1)+(y_t1*y_t1)+(z_t1*z_t1)
    r2_2 = (x_t2*x_t2)+(y_t2*y_t2)+(z_t2*z_t2)
    r5_1 = (r2_1*r2_1) * sqrt(r2_1)
    r5_2 = (r2_2*r2_2) * sqrt(r2_2)
    tmp_1 = 1.0/r5_1
    tmp_2 = 1.0/r5_2
    PCS_1 = (tmp_1*(ax1 *(3.0*z_t1*z_t1 -r2_1)+ rh1*1.5*(x_t1*x_t1- y_t1*y_t1)))
    PCS_2 = (tmp_2*(ax2 *(3.0*z_t2*z_t2 -r2_2)+ rh2*1.5*(x_t2*x_t2- y_t2*y_t2)))
    err = y_v - (PCS_1 + PCS_2)
    return err

def PCSPREMon2SFreeCErr(p, y_v_pcs, y_v_pre, x1 ,y1, z1, wt_pcs, wt_pre):
    """
     FIXME
     @param p:
     @type p:
     @param y_v_pcs:
     @type y_v_pcs:
     @param y_v_pre:
     @type y_v_pre:
     @param x1:
     @type x1:
     @param y1:
     @type y1:
     @param z1:
     @type z1:
     @param wt_pcs:
     @type wt_pcs:
     @param wt_pre:
     @type wt_pre:
    """
    xm1,ym1,zm1, xm2,ym2,zm2, ax1,rh1, ax2,rh2, a1,b1,g1, a2,b2,g2, c = p
    r_v1     = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2     = sqrt((x1-xm2)**2 +(y1-ym2)**2 + (z1-zm2)**2)
    zyz_rot1 = ZYZRot(a1, b1, g1)
    zyz_rot2 = ZYZRot(a2, b2, g2)
    X1, X2   = (x1 - xm1), (x1 - xm2)
    Y1, Y2   = (y1 - ym1), (y1 - ym2)
    Z1, Z2   = (z1 - zm1), (z1 - zm2)
    x_t1     = zyz_rot1[0][0]*X1 + zyz_rot1[0][1]*Y1 + zyz_rot1[0][2]*Z1
    y_t1     = zyz_rot1[1][0]*X1 + zyz_rot1[1][1]*Y1 + zyz_rot1[1][2]*Z1
    z_t1     = zyz_rot1[2][0]*X1 + zyz_rot1[2][1]*Y1 + zyz_rot1[2][2]*Z1
    x_t2     = zyz_rot2[0][0]*X2 + zyz_rot2[0][1]*Y2 + zyz_rot2[0][2]*Z2
    y_t2     = zyz_rot2[1][0]*X2 + zyz_rot2[1][1]*Y2 + zyz_rot2[1][2]*Z2
    z_t2     = zyz_rot2[2][0]*X2 + zyz_rot2[2][1]*Y2 + zyz_rot2[2][2]*Z2
    r2_1     = (x_t1*x_t1)+(y_t1*y_t1)+(z_t1*z_t1)
    r2_2     = (x_t2*x_t2)+(y_t2*y_t2)+(z_t2*z_t2)
    r5_1     = (r2_1*r2_1) * sqrt(r2_1)
    r5_2     = (r2_2*r2_2) * sqrt(r2_2)
    tmp_1    = 1.0/r5_1
    tmp_2    = 1.0/r5_2
    PCS_1    = (tmp_1*(ax1*(3.0*z_t1*z_t1-r2_1)+rh1*1.5*(x_t1*x_t1-y_t1*y_t1)))
    PCS_2    = (tmp_2*(ax2*(3.0*z_t2*z_t2-r2_2)+rh2*1.5*(x_t2*x_t2-y_t2*y_t2)))
    err_pcs  = y_v_pcs - (PCS_1 + PCS_2)
    err_pre  = y_v_pre - (c/r_v1**6 + c/r_v2**6)
    err      = wt_pcs*(err_pcs) + wt_pre*(err_pre)
    return err

def PCSPRE1M2SFC(p0, meas_pcs, meas_pre, x1 ,y1, z1, c, wt_pcs, wt_pre):
    """
     Optimize for two centres using both PCS and PRE data simultaneously
    """
    xm1,ym1,zm1, xm2,ym2,zm2, ax1,rh1, ax2,rh2, a1,b1,g1, a2,b2,g2 = p0
    r_v1     = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2     = sqrt((x1-xm2)**2 +(y1-ym2)**2 + (z1-zm2)**2)
    zyz_rot1 = ZYZRot(a1, b1, g1)
    zyz_rot2 = ZYZRot(a2, b2, g2)
    X1, X2   = (x1 - xm1), (x1 - xm2)
    Y1, Y2   = (y1 - ym1), (y1 - ym2)
    Z1, Z2   = (z1 - zm1), (z1 - zm2)
    x_t1     = zyz_rot1[0][0]*X1 + zyz_rot1[0][1]*Y1 + zyz_rot1[0][2]*Z1
    y_t1     = zyz_rot1[1][0]*X1 + zyz_rot1[1][1]*Y1 + zyz_rot1[1][2]*Z1
    z_t1     = zyz_rot1[2][0]*X1 + zyz_rot1[2][1]*Y1 + zyz_rot1[2][2]*Z1
    x_t2     = zyz_rot2[0][0]*X2 + zyz_rot2[0][1]*Y2 + zyz_rot2[0][2]*Z2
    y_t2     = zyz_rot2[1][0]*X2 + zyz_rot2[1][1]*Y2 + zyz_rot2[1][2]*Z2
    z_t2     = zyz_rot2[2][0]*X2 + zyz_rot2[2][1]*Y2 + zyz_rot2[2][2]*Z2
    r2_1     = (x_t1*x_t1)+(y_t1*y_t1)+(z_t1*z_t1)
    r2_2     = (x_t2*x_t2)+(y_t2*y_t2)+(z_t2*z_t2)
    r5_1     = (r2_1*r2_1) * sqrt(r2_1)
    r5_2     = (r2_2*r2_2) * sqrt(r2_2)
    tmp_1    = 1.0/r5_1
    tmp_2    = 1.0/r5_2
    PCS_1    = (tmp_1*(ax1*(3.0*z_t1*z_t1-r2_1)+rh1*1.5*(x_t1*x_t1-y_t1*y_t1)))
    PCS_2    = (tmp_2*(ax2*(3.0*z_t2*z_t2-r2_2)+rh2*1.5*(x_t2*x_t2-y_t2*y_t2)))
    err_pcs  = meas_pcs - (PCS_1 + PCS_2)
    err_pre  = meas_pre - (c/r_v1**6 + c/r_v2**6)
    err      = wt_pcs*(err_pcs) + wt_pre*(err_pre)
    return err



def InterSpheres(p0, m1,m2,mp,r1,r2,r3):
    """
     Defines the equation for three spheres:
      centered at m1 with radius r1
      centered at m2 with radius r2
      centered at mp with radius r3
     This is used for rootsolving (See: ExplorePara.do_trilat())
    """
    #TODO: Test that his method works
    eq1  = ((p0[0]-m1[0])**2 + (p0[1]-m1[1])**2 + (p0[2]-m1[2])**2) - r1**2
    eq2  = ((p0[0]-m2[0])**2 + (p0[1]-m2[1])**2 + (p0[2]-m2[2])**2) - r2**2
    eq3  = ((p0[0]-mp[0])**2 + (p0[1]-mp[1])**2 + (p0[2]-mp[2])**2) - r3**2
    vals = [eq1,  eq2,  eq3]
    return vals

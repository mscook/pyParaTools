from     numpy import *
from ParaUtils import *

"""Resual error functions evaluated during minimization"""

#NOTE: Function name read as: PRE (1 M)ODEL (1 C)ENTRE (F)IXED (C)
def PRE1M1CFC(p0, meas, c, x,y,z):
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


#NOTE: Function name read as: PRE(1 M)ODEL(1 C)ENTRE (O)PTIMIZE (C)
def PRE1M1SOC(p0, meas, x,y,z):
    """
     Optimize for a single PRE centre and the constant in the PRE equation to
     a single model. Other parameters are as PRE1M1CFC
     @param p0: A numpy array with estimates for the coordinates of the PRE
         centre and an estimate for the constant in the PRE equations
    """
    #NOTE: Believe it agrees with PREfit
    xm,ym,zm,c = p0
    r_v1 = sqrt((x-xm)**2 +(y-ym)**2 + (z-zm)**2)
    err = meas - c/r_v1**6
    return err

#NOTE: Function name read as: PRE(2 M)ODELS(1 C)ENTRE (F)IXED (C)
def PRE2M1SFC(p0, meas, c, x,y,z, x2,y2,z2):
    """
     Optimize for a single PRE centre to 2 models (dimer) with the constant in
     the PRE equation known. Other parameters are as PRE1M1CFC
     @param x2: A numpy array of 2nd x coordinates for the given PRE values
     @param y2: A numpy array of 2nd y coordinates for the given PRE values
     @param z2: A numpy array of 2nd z coordinates for the given PRE values
    """
    xm,ym,zm = p0
    r_v1 = sqrt( (x-xm)**2 +( y-ym)**2 + ( z-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = meas - (c/r_v1**6 + c/r_v2**6)
    return err

#NOTE: Function name read as: PRE(2 M)ODELS(1 C)ENTRE (O)PTIMIZE (C)
def PRE2M1SOC(p0, meas, x,y,z, x2,y2,z2):
    """
     Optimize for a single PRE centre and the constant in the PRE equation to
     2 models.
    """
    #NOTE: Agrees with PREfit
    xm,ym,zm, c = p0
    r_v1 = sqrt(( x-xm)**2 +( y-ym)**2 + ( z-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = meas - (c/r_v1**6 + c/r_v2**6)
    return err

#NOTE: Function name read as: PRE(1 M)ODELS(2 C)ENTRE (F)IXED (C)
def PRE1M2SFC(p0, meas, c, x,y,z):
    """
     Optimize for two PRE centres (dimer case) to a single model with the
     constant in the PRE equation known.
     @param p0:   A numpy array with estimates for the coordinates of each of
        the two PRE centres in the form [<x,y,z>_1  <x,y,z>_2]
    """
    xm, ym, zm, xm2,ym2,zm2 = p0
    r_v1 = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2 = sqrt((x1-xm2)**2 +(y1-ym2)**2 + (z1-zm2)**2)
    err = meas - (c/r_v1**6 + c/r_v2**6)
    return err

def PREMon2SfreecErr(p, y_v, x1,y1,z1):
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
    xm1, ym1, zm1, xm2,ym2,zm2, c = p
    r_v1 = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2 = sqrt((x1-xm2)**2 +(y1-ym2)**2 + (z1-zm2)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PREDimer2SfixedcErr(p, y_v, c, x1,y1,z1, x2,y2,z2):
    """
     FIXME
     @param p:
     @type p:
     @param y_v:
     @type y_v:
     @param c:
     @type c:
     @param x1:
     @type x1:
     @param y1:
     @type y1:
     @param z1:
     @type z1:
     @param x2:
     @type x2:
     @param y2:
     @type y2:
     @param z2:
     @type z2:
    """
    xm1, ym1, zm1, xm2,ym2,zm2 = p
    r_v1 = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2 = sqrt((x2-xm2)**2 +(y2-ym2)**2 + (z2-zm2)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PREDimer2SfreecErr(p, y_v, x1,y1,z1, x2,y2,z2):
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
     @param x2:
     @type x2:
     @param y2:
     @type y2:
     @param z2:
     @type z2:
    """
    xm1, ym1, zm1, xm2,ym2,zm2, c = p
    r_v1 = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2 = sqrt((x2-xm2)**2 +(y2-ym2)**2 + (z2-zm2)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err


loguqmstan1

def PCSMonSSErr(p0, meas, x,y,z):
    """
     FIXME
     @param p0:
     @type p0:
     @param meas:
     @type meas:
     @param x:
     @type x:
     @param y:
     @type y:
     @param z:
     @type z:
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

def PCSPREMon2SFixedCErr(p, y_v_pcs, y_v_pre, x1 ,y1, z1, c, wt_pcs, wt_pre):
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
     @param c:
     @type c:
     @param wt_pcs:
     @type wt_pcs:
     @param wt_pre:
     @type wt_pre:
    """
    xm1,ym1,zm1, xm2,ym2,zm2, ax1,rh1, ax2,rh2, a1,b1,g1, a2,b2,g2 = p
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


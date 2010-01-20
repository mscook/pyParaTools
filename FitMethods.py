from     numpy import *
from ParaUtils import *


def PREMonSSFixedCErr(p, y_v, c, x1,y1,z1):
    xm,ym,zm = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    err = y_v - c/r_v1**6
    return err

def PREMonSSFreeCErr(p, y_v, x1,y1,z1):
    xm,ym,zm,c = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    err = y_v - c/r_v1**6
    return err

def PREDimerSSFixedCErr(p, y_v, c, x1,y1,z1, x2,y2,z2):
    xm,ym,zm = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PREDimerSSFreeCErr(p, y_v, x1,y1,z1, x2,y2,z2):
    xm,ym,zm, c = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PREMon2SFixedCErr(p, y_v, c, x1,y1,z1):
    xm1, ym1, zm1, xm2,ym2,zm2 = p
    r_v1 = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2 = sqrt((x1-xm2)**2 +(y1-ym2)**2 + (z1-zm2)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PREMon2SfreecErr(p, y_v, x1,y1,z1):
    xm1, ym1, zm1, xm2,ym2,zm2, c = p
    r_v1 = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2 = sqrt((x1-xm2)**2 +(y1-ym2)**2 + (z1-zm2)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PREDimer2SfixedcErr(p, y_v, c, x1,y1,z1, x2,y2,z2):
    xm1, ym1, zm1, xm2,ym2,zm2 = p
    r_v1 = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2 = sqrt((x2-xm2)**2 +(y2-ym2)**2 + (z2-zm2)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PREDimer2SfreecErr(p, y_v, x1,y1,z1, x2,y2,z2):
    xm1, ym1, zm1, xm2,ym2,zm2, c = p
    r_v1 = sqrt((x1-xm1)**2 +(y1-ym1)**2 + (z1-zm1)**2)
    r_v2 = sqrt((x2-xm2)**2 +(y2-ym2)**2 + (z2-zm2)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err



def PCSMonSSErr(p, y_v, x1,y1, z1):
    xm, ym, zm, ax, rh, a, b, g = p
    zyz_rot  = ZYZRot(a, b, g)
    X   = x1 - xm
    Y   = y1 - ym
    Z   = z1 - zm
    x_t = zyz_rot[0][0]*X + zyz_rot[0][1]*Y + zyz_rot[0][2]*Z
    y_t = zyz_rot[1][0]*X + zyz_rot[1][1]*Y + zyz_rot[1][2]*Z
    z_t = zyz_rot[2][0]*X + zyz_rot[2][1]*Y + zyz_rot[2][2]*Z
    r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
    r5 = (r2*r2) * sqrt(r2)
    tmp = 1.0/r5
    err = y_v - (tmp*(ax * (3.0*z_t*z_t -r2) + rh*1.5*(x_t*x_t - y_t*y_t)))
    return err

def PCSMon2SErr(p, y_v, x1,y1, z1):
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


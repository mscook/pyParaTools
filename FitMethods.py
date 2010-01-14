from numpy import *
from ParaUtils import *

def PREmfixedcErr(p, y_v, c, x1,y1,z1):
    xm,ym,zm = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    print y_v
    print r_v1
    err = y_v - c/r_v1**6
    return err

def PREmfreecErr(p, y_v, x1,y1,z1):
    xm,ym,zm,c = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    err = y_v - c/r_v1**6
    return err

def PREdfixedcErr(p, y_v, c, x1,y1,z1, x2,y2,z2):
    xm,ym,zm = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PREdfreecErr(p, y_v, x1,y1,z1, x2,y2,z2):
    xm,ym,zm, c = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def PCSmErr(p, y_v, x1,y1, z1):
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


import sys
from numpy import *
import fit_para_helper
import math
##############################################################################
# The core fitting routines
##############################################################################

#####################################################################
# PRE PRE PRE PRE
#####################################################################

###########################################################
# SINGLE METAL SITE OPTIMIZATION
###########################################################

#######################################
# MONOMER
#######################################

def pre_monomer_fixed_c(p, y_v, c, x1,y1,z1):
    # The residual function.
    # Searching for xm, ym, zm
    # c is the constant, p is the current coordinate
    xm,ym,zm = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    err = y_v - c/r_v1**6
    return err

def pre_monomer_free_c(p, y_v, x1,y1,z1):
    # The residual function.
    # Searching for c and xm, ym, zm
    # c is the constant, p is the current coordinate
    xm,ym,zm,c = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    err = y_v - c/r_v1**6
    return err

#######################################
# DIMER
#######################################

def pre_dimer_fixed_c(p, y_v, c, x1,y1,z1, x2,y2,z2):
    # The residual function.
    # Searching for xm, ym, zm
    # c is the constant, p is the current coordinate
    xm,ym,zm = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

def pre_dimer_free_c(p, y_v, x1,y1,z1, x2,y2,z2):
    # The residual function.
    # Searching for xm, ym, zm
    # c is the constant, p is the current coordinate
    xm,ym,zm, c = p
    r_v1 = sqrt((x1-xm)**2 +(y1-ym)**2 + (z1-zm)**2)
    r_v2 = sqrt((x2-xm)**2 +(y2-ym)**2 + (z2-zm)**2)
    err = y_v - (c/r_v1**6 + c/r_v2**6)
    return err

###########################################################
# TWO METAL SITE OPTIMIZATION
###########################################################
#
# TODO
#

#####################################################################
# PCS PCS PCS PCS
#####################################################################

###########################################################
# SINGLE METAL SITE OPTIMIZATION
###########################################################

#######################################
# MONOMER
#######################################

def pcs_monomer(p, y_v, x1,y1,z1):
    ax, rh, xm,ym,zm, a,b,g = p
    PI = math.pi
    scal_const = (1./((12*PI))*10000)
    ax = ax*scal_const
    rh = rh*scal_const
    rot = [["00", "01", "02"], ["10", "11", "12"], ["20", "21", "22"]]
    ca = math.cos(a)
    cb = math.cos(b)
    cg = math.cos(g)
    sa = math.sin(a)
    sb = math.sin(b)
    sg = math.sin(g)
    rot[0][0] = (-sg * sa) + (cb * ca * cg)
    rot[0][1] = ( sg * ca) + (cb * sa * cg)
    rot[0][2] = (-cg * sb)
    rot[1][0] = (-cg * sa) - (cb * ca * sg)
    rot[1][1] = ( cg * ca) - (cb * sa * sg)
    rot[1][2] = ( sg * sb)
    rot[2][0] = ( sb * ca)
    rot[2][1] = ( sb * sa)
    rot[2][2] = ( cb     )
    X = x - xm
    Y = y - ym
    Z = z - zm
    x_t = rot[0][0]*X + rot[0][1]*Y +rot[0][2]*Z
    y_t = rot[1][0]*X + rot[1][1]*Y +rot[1][2]*Z
    z_t = rot[2][0]*X + rot[2][1]*Y +rot[2][2]*Z
    r2  = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
    r5  = (r2*r2) * math.sqrt(r2)
    tmp = 1.0/r5
    pcs = tmp*(ax * (3.0*z_t*z_t -r2) + rh*1.5*(x_t*x_t - y_t*y_t))
    err = yv - pcs
    return err

# Test the ParaData.py object

import sys

sys.path.append("/home/mscook/Desktop/PhD/Projects/pyParaFit")

from ParaData import *

print 80*'-'
print " Testing ParaData.py"
print 80*'-'

################################################################################
#
# Testing the ParaData base class
#
################################################################################

st      = 'H'
sid     = 3
exp_val = 1.11
tol     = 0.1
coord   = [1.546, 2.956, 5.100]
base = ParaData(st, sid, exp_val, tol,coord)


print base.getName()
print base.getName() == st
print base.getId()
print base.getId()   == sid
print base.getVal()
print base.getVal()  == exp_val
print base.getTol()
print base.getTol()  == tol
print base.getCoord()
print base.getCoord() == coord
print base.getCoordx()
print base.getCoordx() ==coord[0]
print base.getCoordy()
print base.getCoordy() ==coord[1]
print base.getCoordz()
print base.getCoordz() ==coord[2]
print base.getType()
print base.getType() =='unknown'
base

st_n      = 'N'
sid_n     = 2
exp_val_n = -0.5
tol_n     = 0.01
coord_n   = [5.000, 11.000, 10.000]
x_n = 5.100
y_n = 11.100
z_n = -10.000

base.setName(st_n)
print base.getName()
base.setId(sid_n)
print base.getId()
base.setVal(exp_val_n)
print base.getVal()
base.setTol(tol_n)
print base.getTol()
base.setCoord(coord_n)
print base.getCoord()
base.setCoordx(x_n)
print base.getCoordx()
base.setCoordy(y_n)
print base.getCoordy()
base.setCoordz(z_n)
print base.getCoordz()
print base.getCoord()


################################################################################
#
# Testing the PCSDataclass
#
################################################################################

st      = 'H'
sid     = 3
exp_val = 1.11
tol     = 0.1
coord   = [1.546, 2.956, 5.100]
pcs = PCSData(st, sid, exp_val, tol,coord)
print pcs.getType()
pcs

################################################################################
#
# Testing the PREDataclass
#
################################################################################

st      = 'H'
sid     = 3
exp_val = 1.11
tol     = 0.1
coord   = [1.546, 2.956, 5.100]
pre = PREData(st, sid, exp_val, tol,coord)
print pre.getType()
pre

################################################################################
#
# Testing the RDCDataclass
#
################################################################################

st      = 'N'
sid     = 3
exp_val = 1.11
tol     = 0.1
coord   = [1.546, 2.956, 5.100]
st2     = 'H'
coord2  = [2.546, 3.956, 6.100]

rdc = RDCData(st, sid, exp_val, tol,coord, st2, coord2)

print rdc.getName()
print rdc.getName() == st+st2
print rdc.getId()
print rdc.getId()   == sid
print rdc.getVal()
print rdc.getVal()  == exp_val
print rdc.getTol()
print rdc.getTol()  == tol
print rdc.getCoord()
a,b= rdc.getCoord()
print a == coord
print b == coord2
print rdc.getCoordx()
a,b = rdc.getCoordx()
print a ==coord[0]
print b == coord2[0]
print rdc.getCoordy()
a,b = rdc.getCoordy()
print a ==coord[1]
print b == coord2[1]
print rdc.getCoordz()
a,b = rdc.getCoordz()
print a ==coord[2]
print b == coord2[2]
print rdc.getType()
print rdc.getType() =='rdc'
a, b = rdc.lookupMGR()
print a,b

#getName
#getCoord(self):
#getCoordx(self):
#getCoordy(self):
#getCoordz(self):
#setName2
#setCoord2(self, coord2):
#setCoord2x(self, coord2_x):
#setCoord2y(self, coord2_y):
#setCoord2z(self, coord2_z):

#        self.s_coord2 = coord2_z
#
#        print rdc.getType()
#rdc.printNice()


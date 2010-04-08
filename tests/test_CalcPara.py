##############################################################################
# Test the CalcPara class
##############################################################################

#TODO: Generate comparative test data using Numbat and PREfit
#TODO: Test all three methods


import sys
sys.path.append("/home/mscook/Desktop/PhD/Projects/pyParaTools")

from   ParaParser  import *
from   CalcPara    import *
from   ExplorePara import *
import os

analysis = ExplorePara()

fail = 0

print 80*'-'
print "Testing CalcPara.py"
print 80*'-'
print "PCS"
pcs_in = ['PROTOCOL_NAME', 'pcs', 'STRUCTURES/epsilon.pdb',
'DATASETS/PCS/EPSILON/PCS_epsilon_CNH.npc', '-5.866', '-0.253', '3.113', '39.994', '4.526', '111.122', '105.981', '110.256']
pcs = PCSParser(pcs_in)
pcs.doParse()
pcs_calcer = CalcPara()
pcs_calcer.PCS(pcs, 'ZYZ')
analysis.buildNumbatTBL(pcs, 'test_pcs.npc')
print "Comparing with known..."
exact  = open('DATASETS/PCS/EPSILON/EPS_TESTING.npc').readlines()
calced = open('test_pcs.npc').readlines()
if len(exact) != len(calced):
    print 'FAILED'
    fail+=1
for i in range(0, len(exact)):
    a = exact[i].split()
    b = calced[i].split()
    # Tol set to deal with rounding in Numbat
    if float(a[2]) - float(b[2]) > 0.0011:
        print 'FAILED'
        print b[0], float(a[2]) - float(b[2])
        fail+=1
print "Number of failures", fail
print 80*'-'
fail = 0


print "RDC"
rdc_in = ['PROTOCOL_NAME', 'rdc', 'STRUCTURES/epsilon.pdb',
'DATASETS/RDC/epsilon.rdc', '-23.392', '37.731', '57.846', '25.591', '230.403']
rdc = RDCParser(rdc_in)
rdc.doParse()
rdc_calcer = CalcPara()
rdc_calcer.RDC(rdc, 'ZYZ')
analysis.buildNumbatTBL(rdc, 'test_rdc.npc')
print "Comparing with known..."
exact  = open('DATASETS/RDC/EPS_TESTING.npc').readlines()
calced = open('test_rdc.npc').readlines()
if len(exact) != len(calced):
    print 'FAILED'
    fail+=1
for i in range(0, len(exact)):
    a = exact[i].split()
    b = calced[i].split()
    # Tolerance here is 0.03
    # This is to do with slight differences in constants
    if float(a[2]) - float(b[2]) > 0.03:
        print 'FAILED'
        print b[0], float(a[2]) - float(b[2])
        fail+=1
print "Number of failures", fail
print 80*'-'
fail = 0


print "PRE"
pre_in = ['PROTOCOL_NAME', 'pre', 'STRUCTURES/epsilon.pdb',
'DATASETS/PRE/EPS_TESTING_PRE.npc', '-5.866', '-0.253', '3.113', '680000000']
pre = PREParser(pre_in)
pre.doParse()
pre_calcer = CalcPara()
pre_calcer.PRE(pre)
analysis.buildNumbatTBL(pre, 'test_pre.npc')
print "Comparing with known..."
exact  = open('DATASETS/PRE/EPS_TESTING_PRE.npc').readlines()
calced = open('test_pre.npc').readlines()
if len(exact) != len(calced):
    print 'FAILED'
    fail+=1
for i in range(0, len(exact)):
    a = exact[i].split()
    b = calced[i].split()
    if float(a[2]) - float(b[2]) >= 0.005:
        print 'FAILED'
        print b[0], float(a[2]) - float(b[2])
        fail+=1
print "Number of failures", fail
print 80*'-'
fail = 0

os.system("rm *.npc")


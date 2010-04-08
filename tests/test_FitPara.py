#Test the FitPara classes

import sys

sys.path.append("/home/mscook/Desktop/PhD/Projects/pyParaTools")

from ParaParser import *
from CalcPara import *
from FitPara import *
from ParaUtils import *

print 80*'-'
print " Testing FitPara.py"
print 80*'-'


#python test_FitPara.py pre ~/Desktop/PREfit_python/TESTDATA/PDB/m0.pdb ~/Desktop/PREfit_python/TESTDATA/PRE/EXPERIMENTAL/PREdata_intra_contribution.pre -4.20536  36.52032 -3.35519 13554627736.87463
#pre1 = PREParser(sys.argv)
#pre1.doParse()
#pre_calc = CalcPara()
#pre_calc.PRE(pre1)
#pre_fit  = FitPara()
#pre_fit.pre_monomer_fixed_c(pre1)

##python test_FitPara.py pcs STRUCTURES/epsilon.pdb DATASETS/PCS/PCS_epsilon_CNH.npc 1 1 1 30 10 90 90 90
#pcs1 = PCSParser(sys.argv)
#pcs1.doParse()
##print pcs1.getParsed()
#pcs_calcer = CalcPara()
#pcs_calcer.PCSZYZ(pcs1)
#print pcs_1
#pcs_fit  = FitPara()
#pcs_fit.pcs_monomer(pcs1)


rdc_in = ['PROTOCOL_NAME', 'rdc', 'STRUCTURES/epsilon.pdb',
'DATASETS/RDC/epsilon.rdc','-25.392', '39.731', '59.846', '30.591', '240.403']
print '-23.392', '37.731', '57.846', '25.591', '230.403'
rdc = RDCParser(rdc_in)
rdc.doParse()
#rdc_calcer = CalcPara()
#rdc_calcer.RDC(rdc, 'ZYZ')
B0    = rdc.getB0()
temp  = rdc.getTemp()
S     = rdc.getOrder()
g1    = lookupMGR('H')
g2    = lookupMGR('N')[1]
scal  = rdcScal(S, g1, g2, B0, temp)

rdc_fitter = FitPara()
rdc_fitter.RDC(rdc, 0, scal)


from ParaParser  import *
from CalcPara    import *
from FitPara     import *
from ExplorePara import *
from ParaUtils   import *

print
print
print """     pyParaFit - fitting paramagnetic observables to theoretical models
           pyParaFit  Copyright (C) 2010  Mitchell J Stanton-Cook """
print 79*'-'
print """pyParaFit comes with ABSOLUTELY NO WARRANTY; for details read LICENCE
This is free software, and you are welcome to redistribute it
under certain conditions; read LICENCE for details. """
print 79*'-'
print


################################################################################
# EXAMPLES
################################################################################

in_1 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
'tests/DATASETS/PCS/1.npc', '-3.', '-2.', '9.','-15.0', '15.0', '140', '15', '81']
#in_1 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/1.npc', '-5.82', '-4.10', '7.30','-16.98', '17.41',
#'139.09', '14.36', '81.85']


#pcs_1 = PCSParser(in_1)
#pcs_1.doParse()
#pcs1_calcer = CalcPara()
#pcs1_calcer.PCS(pcs_1, 'ZYZ')
#analysis = ExplorePara()
#analysis.paraDevs(pcs_1)
#pcs_fit  = FitPara()
#res = pcs_fit.PCS(pcs_1, 0)
#print "x,y,z:", res[0], res[1], res[2]
#print "Ax, Rh", FromVVU(res[3]), FromVVU(res[4])
#print res[5], res[6], res[7]
#print "A,B,G",  FixAngle(res[5]), FixAngle(res[6]), FixAngle(res[7])
##------------
pre_in = ['pyParaFit.py', 'pre', 'tests/STRUCTURES/2ZTAa_H_only.pdb',
'tests/DATASETS/PRE/PREdata_intra_contribution.pre', '0', '0', '0', '0']
pre_1_m1 = PREParser(pre_in)
pre_1_m1.doParse()

pre_1_m2 = PREParser(pre_in)
pre_1_m2.setModel(1)
pre_1_m2.doParse()

pre_1_calcer = CalcPara()
pre_1_calcer.PRE(pre_1_m1)

analysis = ExplorePara()
analysis.paraSummary(pre_1_m1)
#analysis.plotParaDevs(parsed)
pre_fit_calcer = FitPara()
#res = pre_fit_calcer.PRE(pre_1_m1, 3, pre_1_m2)
res = pre_fit_calcer.PRE(pre_1_m1, 1)
print "x,y,z,c", res[0], res[1], res[2], res[3]
opt_metal = zeros(3)
opt_metal[0], opt_metal[1], opt_metal[2] = res[0], res[1], res[2]
opt_c = zeros(1)
opt_c[0] = res[3]
pre_1_m1.setMetalLoc(opt_metal)
pre_1_m1.setConstant(opt_c)
pre_1_calcer.PRE(pre_1_m1)
analysis.paraSummary(pre_1_m1)




#in_2 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/PCS_epsilon_CNH.npc','1','1','1', '30','10','90','90','90']

#in_3 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/PCS_epsilon_CNH.npc','1','1','1', '30','10','90','90','90']

#in_4 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/PCS_epsilon_CNH.npc','1','1','1', '30','10','90','90','90']


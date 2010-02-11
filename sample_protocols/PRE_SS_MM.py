import sys

sys.path.append("/home/mscook/Desktop/PhD/Projects/pyParaFit")

from ParaParser  import *
from CalcPara    import *
from FitPara     import *
from ExplorePara import *
from ParaUtils   import *


################################################################################
# SAMPLE PROTOCOL: PRE SingleSite, MultipleModels
#
# Fit the location of the PRE centre using a symmetric dimer structure.
# Look at both fixed and free PRE constant fits. Does some simple analysis
#
# 01/10 Mitchell Stanton-Cook
################################################################################
print
print
# Setup for the type of data, the structural model, dataset and initial search
# parameters
pre_in = ['pyParaFit.py', 'pre', '../tests/STRUCTURES/2ZTAa_H_only.pdb',
'../tests/DATASETS/PRE/PREdata_intra_contribution.pre', '-7', '-26', '3', '3500000000']

#Build a parser object and parse
pre_1_m1 = PREParser(pre_in)
pre_1_m1.doParse()
# Now lets setup for the second model
pre_1_m2 = PREParser(pre_in)
pre_1_m2.setModel(1)
pre_1_m2.doParse()
#Now lets perform a non-linear fit to the "monomer" to determine the position
#of the metal ion
pre_fitter = FitPara()
# The additional argument '1' in this case specifies that we want to use a
# single model and optimize c
res = pre_fitter.PRE(pre_1_m1, 3, pre_1_m2)
print "Optimization with initial params:", pre_1_m1.getMetalLoc(), pre_1_m1.getConstant()
print
print "Determined x,y,z and c", res[0], res[1], res[2], res[3]
print
# Now lets update the metal and c value to the optimized values
opt_metal = zeros(3)
opt_metal[0], opt_metal[1], opt_metal[2] = res[0], res[1], res[2]
opt_c = zeros(1)
opt_c[0] = res[3]
pre_1_m1.setMetalLoc(opt_metal)
pre_1_m1.setConstant(opt_c)
pre_1_m2.setMetalLoc(opt_metal)
pre_1_m2.setConstant(opt_c)
#  Now re-calculate
pre_calcer = CalcPara()
pre_calcer.PRE(pre_1_m1)
pre_calcer.PRE(pre_1_m2)
# Do some analysis
print "Analysis of PRE agreement for fit to dimer with c optimized"
analysis = ExplorePara()
analysis.paraSummaryMulti(pre_1_m1, pre_1_m2)

# Now, lets look at an optimization with fixed c
init_metal = zeros(3)
fixed_c = zeros(1)
init_metal[0], init_metal[1], init_metal[2] = -7.000, -26.000, 3.000
fixed_c[0] = 680000000
pre_1_m1.setMetalLoc(init_metal)
pre_1_m1.setConstant(fixed_c)
pre_1_m2.setMetalLoc(init_metal)
pre_1_m2.setConstant(fixed_c)
res = pre_fitter.PRE(pre_1_m1, 2, pre_1_m2)
print "Optimization with c fixed to", fixed_c
print "Determined x,y,z and c", res[0], res[1], res[2]
opt_metal = zeros(3)
opt_metal[0], opt_metal[1], opt_metal[2] = res[0], res[1], res[2]
pre_1_m1.setMetalLoc(opt_metal)
pre_1_m1.setConstant(fixed_c)
pre_1_m2.setMetalLoc(opt_metal)
pre_1_m2.setConstant(fixed_c)
#  Now re-calculate
pre_calcer = CalcPara()
pre_calcer.PRE(pre_1_m1)
pre_calcer.PRE(pre_1_m2)
analysis.paraSummaryMulti(pre_1_m1, pre_1_m2)

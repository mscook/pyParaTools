import sys

sys.path.append("/home/mscook/Desktop/PhD/Projects/pyParaFit")

from ParaParser  import *
from CalcPara    import *
from FitPara     import *
from ExplorePara import *
from ParaUtils   import *


################################################################################
# SAMPLE PROTOCOL: PRE SingleSite, SingleModel
#
# Fit the location of the PRE centre using a single structure. Look at both
# fixed and free PRE constant fits. Does some simple analysis
#
# 01/10 Mitchell Stanton-Cook
################################################################################

# Setup for the type of data, the structural model, dataset and initial search
# parameters
pre_in = ['pyParaFit.py', 'pre', '../tests/STRUCTURES/2ZTAa_H_only.pdb',
'../tests/DATASETS/PRE/PREdata_intra_contribution.pre', '0', '0', '0', '0']

#Build a parser object and parse
pre_1 = PREParser(pre_in)
pre_1.doParse()
#Lets calculate the PRE effect (by building a calcer object)experienced given
#these initial search params
pre_1_calcer = CalcPara()
pre_1_calcer.PRE(pre_1)
#Calcing simply populates the data structure with calculated values. We really
#need to explore (analyze our dataset))
analysis = ExplorePara()
analysis.paraSummary(pre_1)
#analysis.paraSummaryPlot(pre_1)
#As we can see our initial search params are poor. Now lets perform a non-linear
#fit to find more "optimum values""
pre_1_fitter = FitPara()
# The additional argument '1' in this case specifies that we want to use a
# single model and optimize c
res = pre_1_fitter.PRE(pre_1, 1)
print "Determined x,y,z and c", res[0], res[1], res[2], res[3]
# Now lets update the metal and c value to the optimized values
opt_metal = zeros(3)
opt_metal[0], opt_metal[1], opt_metal[2] = res[0], res[1], res[2]
opt_c = zeros(1)
opt_c[0] = res[3]
pre_1.setMetalLoc(opt_metal)
pre_1.setConstant(opt_c)
# Now re-calculate
pre_1_calcer.PRE(pre_1)
# Get the summary and the Qfactor (quality of the agreement of obs with clac)
analysis.paraSummary(pre_1)
analysis.paraSummaryPlot(pre_1)
print "Qfactor is:", analysis.QFactor(pre_1)

#Now, just check the agreement using optimization with fixed c
opt_metal = zeros(3)
opt_metal[0], opt_metal[1], opt_metal[2] = 0.0, 0.0, 0.0
pre_1.setMetalLoc(opt_metal)
pre_1_calcer.PRE(pre_1)
# Get the summary and the Qfactor (quality of the agreement of obs with clac)
analysis.paraSummary(pre_1)
res = pre_1_fitter.PRE(pre_1, 0)
print "Determined x,y,z ", res[0], res[1], res[2]


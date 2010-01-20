from ParaParser import *
from CalcPara   import *
from FitPara    import *

print 80*'-'
print " pyParaFit.py"
print 80*'-'

# Get the data
# For datasets simplest to pass as arguments
#in_data = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/PCS_epsilon_CNH.npc','1','1','1', '30','10','90','90','90']

#Build a parser and parse
#pcs_eps = PCSParser(in_data)
#pcs_eps.doParse()

# Build a calcer and calc PCS
#pcs_calcer = CalcPara()
#pcs_calcer.PCSZYZ(pcs_eps)

################################################################################
# EXAMPLE
################################################################################

in_1 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
'tests/DATASETS/PCS/1.npc', '-5.', '-4.', '7.','-15.0', '15.0',
'140', '15', '81']
#in_1 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/1.npc', '-5.82', '-4.10', '7.30','-16.98', '17.41',
#'139.09', '14.36', '81.85']

#Build a parser & parse
pcs_1 = PCSParser(in_1)
pcs_1.doParse()
#Build a calcer & calc
#pcs1_calcer = CalcPara()
#pcs1_calcer.PCS(pcs_1, 'ZYZ')
#Get the Obs, Calced
vals = pcs_1.getParsed()
a,b,c,d,e = 'Atom Type', '   Residue Number', 'Experimental', 'Calculated', 'Deviation'
print '%s%12s%15s%15s%15s' % (a,b,c,d,e)
for pObject in range(0, len(vals)):
    print '%s%20i%20.3f%20.3f%20.3f' % (vals[pObject].getName(), vals[pObject].getId(), vals[pObject].getVal(), vals[pObject].getCVal(), vals[pObject].getVal()-vals[pObject].getCVal())

pcs_fit_1  = FitPara()
pcs_fit_1.pcs_monomer(pcs_1)


#in_2 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/PCS_epsilon_CNH.npc','1','1','1', '30','10','90','90','90']

#in_3 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/PCS_epsilon_CNH.npc','1','1','1', '30','10','90','90','90']

#in_4 = ['pyParaFit.py', 'pcs', 'tests/STRUCTURES/epsilon.pdb',
#'tests/DATASETS/PCS/PCS_epsilon_CNH.npc','1','1','1', '30','10','90','90','90']


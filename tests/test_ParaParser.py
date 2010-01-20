import sys


sys.path.append("/home/mscook/Desktop/PhD/Projects/pyParaFit")

from ParaParser import *

print 80*'-'
print " Testing ParaParser.py"
print 80*'-'

#python test_ParaParser.py pcs ~/Desktop/PREfit_python/TESTDATA/PDB/m0.pdb ~/Desktop/PREfit_python/TESTDATA/PCS/EXPERIMENTAL/LeucineZipper_3.0YbDPA3.npc 1 1 1 30 10 90 90 90
###pcs1 = PCSParser(sys.argv)
###print  pcs1.getNumModels()
###pcs1.doParse()
###print pcs1.getParsed()

#

pre1 = PREParser(sys.argv)
pre1.doParse()
print pre1.getParsed()
#
#rdc = RDCParser(sys.argv)
#rdc.doParse()
#print rdc.parsed


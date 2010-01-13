import sys


sys.path.append("/home/mscook/Desktop/PhD/Projects/pyParaFit")

from ParaParser import *

print 80*'-'
print " Testing ParaParser.py"
print 80*'-'


parse = ParaParser(sys.argv)
a = parse.getNumModels()
print a
parse.doParse()
print parse.getParsed()

rdc = RDCParser(sys.argv)
rdc.doParse()
print rdc.parsed


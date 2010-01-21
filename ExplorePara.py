"""A class to explore/visualize paramagnetic datasets"""

class ExplorePara:

    def paraDevs(self ,ParsedObj):
        data     = ParsedObj.getParsed()
        print 'AtomType  ResidueNumber  Experimental  Calculated  Deviation'
        for pObject in range(0, len(data)):
            print '%s%20i%20.3f%20.3f%20.3f' % (data[pObject].getName(), data[pObject].getId(), data[pObject].getVal(), data[pObject].getCVal(), data[pObject].getVal()-data[pObject].getCVal())

    def plotParaDevs(self, ParsedObj):
        from numpy import *
        import pylab
        data     = ParsedObj.getParsed()
        obsv, calc = [], []
        for pObject in range(0, len(data)):
            obsv.append(data[pObject].getVal())
            calc.append(data[pObject].getCVal())
        pylab.plot(array(obsv), array(calc))
        pylab.xlabel('Observed PCS [ppm]')
        pylab.ylabel('Calced PCS [ppm]')
        pylab.title('PCS agreement')
        pylab.show()


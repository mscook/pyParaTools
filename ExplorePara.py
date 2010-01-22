from numpy import *

"""A class to explore/visualize paramagnetic datasets"""

class ExplorePara:
#TODO: Add numerical spins
#TODO: Add SFP
#TODO: Add RDC params
#TODO: Add Qfactor/R-free
#TODO: Add isosurface

    def paraSummary(self ,ParsedObj):
        data     = ParsedObj.getParsed()
        print
        print 80*'-'
        print 'ATOM  RES  EXP         CALC1        DEV'
        print 80*'-'
        for pObject in range(0, len(data)):
            aname = data[pObject].getName()
            rnum  = data[pObject].getId()
            exp_v = data[pObject].getVal()
            cal_v = data[pObject].getCVal()
            dev   = data[pObject].getVal()-data[pObject].getCVal()
            print '%s%6i%12.3f%12.3f%12.3f' % (aname,rnum,exp_v,cal_v,dev)
        print 80*'-'

    def paraSummaryMulti(self ,ParsedObj1, ParsedObj2):
        data1     = ParsedObj1.getParsed()
        data2     = ParsedObj2.getParsed()
        print
        print 80*'-'
        out1 = 'ATOM  RES  EXP         CALC1     CALC2      CALCT        DEV'
        out2 = '              %S2'
        print out1+out2
        print 80*'-'
        for pObject in range(0, len(data1)):
            aname = data1[pObject].getName()
            rnum  = data1[pObject].getId()
            exp_v = data1[pObject].getVal()
            cal_1 = data1[pObject].getCVal()
            cal_2 = data2[pObject].getCVal()
            cal_t = data1[pObject].getCVal()+data2[pObject].getCVal()
            dev   = cal_t - exp_v
            pc_s2 = cal_2/cal_t
            print '%s%6i%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f' % \
            (aname, rnum, exp_v, cal_1, cal_2, cal_t, dev, pc_s2)
        print 80*'-'


    def paraSummaryPlot(self, ParsedObj):
        import pylab
        data     = ParsedObj.getParsed()
        obsv, calc = [], []
        for pObject in range(0, len(data)):
            obsv.append(data[pObject].getVal())
            calc.append(data[pObject].getCVal())
        pylab.plot(array(obsv), array(calc))
        pylab.xlabel('Observed ')
        pylab.ylabel('Calced')
        pylab.title('Fit agreement')
        pylab.show()

    def QFactor(self, ParsedObj):
        import math
        sum1 = 0.0
        sum2 = 0.0
        data     = ParsedObj.getParsed()
        for pObject in range(0, len(data)):
            obs  = data[pObject].getVal()
            calc = data[pObject].getCVal()
            sum1 = sum1 + (obs - calc)**2
            sum2 = sum2 + (obs)**2
        Qfac =  math.sqrt(sum1/sum2)
        return Qfac


from numpy import *
from scipy import *
import math

"""A class to explore/visualize paramagnetic datasets"""

class ExplorePara:
#TODO: Add numerical spins
#TODO: Add SFP
#TODO: Add RDC params
#TODO: Add isosurface

    def paraSummary(self ,ParsedObj):
        data     = ParsedObj.getParsed()
        print
        print 80*'-'
        print 'ATOM  RES  EXP         CALC1        DEV           DEV/EXP'
        print 80*'-'
        for pObject in range(0, len(data)):
            aname  = data[pObject].getName()
            rnum   = data[pObject].getId()
            exp_v  = data[pObject].getVal()
            cal_v  = data[pObject].getCVal()
            dev    = data[pObject].getVal()-data[pObject].getCVal()
	    perDev = abs(dev)/exp_v	
            print '%s%6i%12.3f%12.3f%12.3f%12.3f' % (aname,rnum,exp_v,cal_v,dev, perDev)
        print 80*'-'

    def paraSummaryMulti(self ,ParsedObj1, ParsedObj2):
    #NOTE: Check ie. -ve %'s and > 1
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


    def paraSummaryPlotMulti(self, ParsedObj1, ParsedObj2):
        import pylab
        data1     = ParsedObj1.getParsed()
        data2     = ParsedObj2.getParsed()
        obsv, calc = [], []
        for pObject in range(0, len(data1)):
            obsv.append(data[pObject].getVal())
            calc_t = data1[pObject].getCVal()+data2[pObject].getCVal()
            calc.append(calc_t)
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


    def QFactorMulti(self, ParsedObj1, ParsedObj2):
        import math
        sum1 = 0.0
        sum2 = 0.0
        data1     = ParsedObj1.getParsed()
        data2     = ParsedObj2.getParsed()
        for pObject in range(0, len(data1)):
            obs  = data1[pObject].getVal()
            calc_t = data1[pObject].getCVal()+data2[pObject].getCVal()
            sum1 = sum1 + (obs - calc_t)**2
            sum2 = sum2 + (obs)**2
        Qfac =  math.sqrt(sum1/sum2)
        return Qfac



    def describeData(self, ParsedObj):
        exp = ParsedObj.getAllMeasarray()
        print mean(exp)
        print var(exp)


    def __buildPymolScript(self, curfname, m1m2vec, mpvec):
        f_out    = open('do_trilat.pml', 'w')
        f        = curfname.split('/')
        f_clean  = f[len(f)-1][:-4]
        target_n = f_clean+"_soln"
        f_out.write("copy "+target_n+", "+f_clean+"\n")
        f_out.write("pair_fit \\"+"\n")
        f_out.write(f_clean+"_soln"+" & i. 500, ALIGN_SOLN & i. 501, \\"+"\n")
        f_out.write(f_clean+"_soln"+" & i. 901, ALIGN_SOLN & i. 900, \\"+"\n")
        f_out.write(f_clean+"_soln"+" & i. 501, ALIGN_SOLN & i. 500, "+"\n")
        i = 0
        f_out.write('save 0.pdb, '+target_n+'\n')
        while i <= 355:
            i = i+5
            f_out.write("copy "+target_n+"_"+str(i) +", "+target_n+"\n")
            f_out.write("rotate ["+str(m1m2vec[0])+","+str(m1m2vec[1])+"," \
            +str(m1m2vec[2])+"], "+str(i)+", "+target_n+"_"+str(i) +        \
            ", camera=0, origin=["+str(mpvec[0])+","+str(mpvec[1])+","+  \
            str(mpvec[2])+"]"+"\n")
            f_out.write('save '+str(i)+".pdb, "+target_n+"_"+str(i)+'\n')
        f_out.close()


    def __trilatPDB(self, m1,m2,sk,soln):
	#TODO: Need to add the metal sites to the structure file.
        f_out    = open('ALIGN_SOLN.pdb', 'w')
        a, aid, at, rt, mod, resn = \
          "ATOM", [5000,5001,9000, 9001], ["M1","M2","SU","SK"], \
          ["MET","UKS","KNS"], "A", [500,501,900,901]
        f_out.write('%s%7i%4s%5s%2s%4i%12.3f%8.3f%8.3f' \
          %(a,aid[0],at[0],rt[0],mod,resn[0], m1[0], m1[1], m1[2])+"\n")
        f_out.write('%s%7i%4s%5s%2s%4i%12.3f%8.3f%8.3f' \
          %(a,aid[1],at[1],rt[0],mod,resn[1], m2[0], m2[1], m2[2])+"\n")
        print "ADD THE FOLLOWING ATOM TO THE APPROPRIATE STRUCTURE/MODEL:"
        print '%s%7i%4s%5s%2s%4i%12.3f%8.3f%8.3f' \
          %(a,aid[3],at[3],rt[2],mod,resn[3], sk[0], sk[1], sk[2])
        f_out.write('%s%7i%4s%5s%2s%4i%12.3f%8.3f%8.3f' \
          %(a, aid[2],at[2],rt[1],mod,resn[2],soln[0],soln[1],soln[2]) +"\n")
        print "THEN IN pyMOL RUN do_trilat.pml"


    def buildPossibleComplex(self, ParsedObj1, ParsedObj2):
        #TODO: Test that his method works
        from scipy.optimize import fsolve
        from FitMethods import InterSpheres
        import os

        m1     = ParsedObj1.getMetalLoc()
        m2     = ParsedObj2.getMetalLoc()
        #Always works with the 1st assigned coordinate
        sk = ParsedObj1.getParsed()[0].getCoord()
        # Do some math
        mpx    = (m1[0]+m2[0])/2.0
        mpy    = (m1[1]+m2[1])/2.0
        mpz    = (m1[2]+m2[2])/2.0
        mp     = zeros(3)
        mp[0], mp[1], mp[2] =  mpx,  mpy, mpz
        m1m2   = zeros(3)
        m1m2   = m1-m2
        # Calculate distances
        dm1_s1 = math.sqrt((m1[0]-sk[0])**2 +(m1[1]-sk[1])**2 +(m1[2]-sk[2])**2)
        dm2_s1 = math.sqrt((m2[0]-sk[0])**2 +(m2[1]-sk[1])**2 +(m2[2]-sk[2])**2)
        dm1_m2 = math.sqrt((m1[0]-m2[0])**2 +(m1[1]-m2[1])**2 +(m1[2]-m2[2])**2)
        dmp_s1 = math.sqrt((mp[0]-sk[0])**2 +(mp[1]-sk[1])**2 +(mp[2]-sk[2])**2)
        # Distances for r2/r1 are switched.
        r2 = dm1_s1
        r1 = dm2_s1
        r3 = dmp_s1
        p0 = \
         [random.randint(-50,50),random.randint(-50,50), random.randint(-50,50)]
        soln = fsolve(InterSpheres, p0, args=(m1,m2,mp,r1,r2,r3), full_output=1)
        check = soln[3]
        if check[0:23] == 'The solution converged.':
            self.__trilatPDB(m1,m2,sk, soln[0])
            self.__buildPymolScript(ParsedObj1.getPDBFn(), m1m2, mp)
            os.system("cp "+ParsedObj1.getPDBFn()+" .")
        else:
            print "Failed. Retry with different initial conditions"


from numpy          import *
from scipy.optimize import leastsq
from FitMethods     import *
from ParaUtils      import *

class FitPara:
#    def __init__(self, pds1, pds2=[], pds3=[], pds4=[], pds5=[], pds6=[],
#     pds7=[], pds8=[], pds9=[], pds10=[], pds11=[], pds12=[]):

#        self.pD1 = pds1.getParsed()
#        xlist1, ylist1, zlist1, expvals_list1 = [], [], [], []
#        for pObject in range(0, len(self.pD)):
#            x,y,z, meas = data[pObject].getCoordx(), data[pObject].getCoordy(), data[pObject].getCoordz(), data[pObject].getVal()
#            xl.append(x)
#            yl.append(y)
#            zl.append(z)
#            measl.append(meas)
#        xl = array(xl)
#        yl = array(yl)
#        zl = array(zl)
#        measl = array(measl)


    def PRE(self, parsed_pre_data, type_of_fit):
        data = parsed_pre_data.getParsed()
        xl, yl, zl, measl = [], [], [], []
        for pObject in range(0, len(data)):
            x,y,z, meas = data[pObject].getCoordx(), data[pObject].getCoordy(), data[pObject].getCoordz(), data[pObject].getVal()
            xl.append(x)
            yl.append(y)
            zl.append(z)
            measl.append(meas)
        xl = array(xl)
        yl = array(yl)
        zl = array(zl)
        measl = array(measl)
        if type_of_fit == 0:
            p0   = parsed_pre_data.getMetalLoc()
            c    = parsed_pre_data.getConstant()
            soln,cov,info,mesg,success = leastsq(PREmfixedcErr, p0, args=(measl, c,xl,yl,zl), full_output=1)
            print soln
            #fit_para_helper.get_fit_summary(soln,cov,info,mesg,success, meas, p0, name)

    def pcs_monomer(self, parsed_pcs_data):
        p0_tmp   = (parsed_pcs_data.getMetalLoc())
        p0 = p0_tmp.tolist()
        p0.append(parsed_pcs_data.getAxialVVU())
        p0.append(parsed_pcs_data.getRhombicVVU())
        p0.append(parsed_pcs_data.getAlpha())
        p0.append(parsed_pcs_data.getBeta())
        p0.append(parsed_pcs_data.getGamma())
        data = parsed_pcs_data.getParsed()
        print p0
        xl, yl, zl, measl = [], [], [], []
        for pObject in range(0, len(data)):
            x,y,z, meas = data[pObject].getCoordx(), data[pObject].getCoordy(), data[pObject].getCoordz(), data[pObject].getVal()
            xl.append(x)
            yl.append(y)
            zl.append(z)
            measl.append(meas)
        xl = array(xl)
        yl = array(yl)
        zl = array(zl)
        measl = array(measl)
        soln,cov,info,mesg,success = leastsq(PCSMonSSErr, p0, args=(measl,xl,yl,zl), full_output=1)
        print soln,cov,info,mesg,success
        #print 'Metal: ', soln[0], soln[1], soln[2]
        sAx = FromVVU(soln[3])
        sRh = FromVVU(soln[4])
        #print 'Ax, Rh:', sAx, sRh
        #print 'Alpha, Beta, Gamma:', soln[5], soln[6], soln[7]
        #fit_para_helper.get_fit_summary(soln,cov,info,mesg,success, meas, p0, name)


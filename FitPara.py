from numpy          import *
from scipy.optimize import leastsq
from FitMethods     import *
from ParaUtils      import *

class FitPara:


    def PCS(self, parsedObj, type_of_fit):
        #TODO: Optimize me please
        x    = parsedObj.getAllXarray()
        y    = parsedObj.getAllYarray()
        z    = parsedObj.getAllZarray()
        meas = parsedObj.getAllMeasarray()
        met = parsedObj.getMetalLoc()
        ax  = parsedObj.getAxialVVU()
        rh  = parsedObj.getRhombicVVU()
        a   = parsedObj.getAxialVVU()
        b   = parsedObj.getAxialVVU()
        g   = parsedObj.getAxialVVU()
        if type_of_fit == 0:
            p0 = zeros(8)
            p0[0], p0[1], p0[2], p0[3], p0[4], p0[5], p0[6], p0[7] = met[0], met[1], met[2], ax, rh, a, b ,g
            soln,cov,info,mesg,success = leastsq(PCSMonSSErr, p0, args=(meas, x,y,z), full_output=1)
            return soln

    #TODO: Optimize me please
    def PRE(self, parsedObj, type_of_fit, parsedObj2="None" ):
        #NOTE: Agrees with PREfit (by fixing c to free optimized value.
        x    = parsedObj.getAllXarray()
        y    = parsedObj.getAllYarray()
        z    = parsedObj.getAllZarray()
        meas = parsedObj.getAllMeasarray()
        c    = parsedObj.getConstant()
        if type_of_fit == 0:
            p0   = parsedObj.getMetalLoc()
            soln,cov,info,mesg,success = leastsq(PRE1M1SFC, p0, args=(meas, c,x,y,z), full_output=1)
        if type_of_fit == 1:
            #NOTE: Agrees with PREfit
            met   = parsedObj.getMetalLoc()
            p0 = zeros(4)
            p0[0], p0[1], p0[2], p0[3] = met[0],met[1],met[2], c[0]
            soln,cov,info,mesg,success = leastsq(PRE1M1SOC, p0, args=(meas, x,y,z), full_output=1)
        if type_of_fit == 2:
            #NOTE: Agrees with PREfit
            x2    = parsedObj2.getAllXarray()
            y2    = parsedObj2.getAllYarray()
            z2    = parsedObj2.getAllZarray()
            p0    = parsedObj2.getMetalLoc()
            soln,cov,info,mesg,success = leastsq(PRE2M1SFC, p0, args=(meas, c,x,y,z,x2,y2,z2), full_output=1)
        if type_of_fit == 3:
            #NOTE: Agrees with PREfit
            x2    = parsedObj2.getAllXarray()
            y2    = parsedObj2.getAllYarray()
            z2    = parsedObj2.getAllZarray()
            met   = parsedObj.getMetalLoc()
            p0 = zeros(4)
            p0[0], p0[1], p0[2], p0[3] = met[0],met[1],met[2], c[0]
            soln,cov,info,mesg,success = leastsq(PRE2M1SOC, p0, args=(meas, x,y,z,x2,y2,z2), full_output=1)

        if success == 1:
            return soln
        else:
            msg1 = "Optimization failed. Please try again with different "
            msg2 = "different intial starting guesses"
            print msg1+msg2


    def RDC(self, parsedObj, type_of_fit):
        # NOTE: This method will come in future- when RDC calc is known to work
        pass


    def PCS_PRE(self, parsedObj, type_of_fit):
        #TODO: Write this method
        pass


    def PCS_RDC(self, parsedObj, type_of_fit):
        # NOTE: This method will come in future (when RDC calc is known to work
        pass


    def PRE_RDC(self, parsedObj, type_of_fit):
        # NOTE: This method will come in future (when RDC calc is known to work
        pass


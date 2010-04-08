from numpy          import *
from scipy.optimize import leastsq
from FitMethods     import *
from ParaUtils      import *

class FitPara:


    def PCSOptions(self):
        """
        Prints information on the types of fits that can be achieved from PCS
        data
        """
        print "0       1 X-t        1 Model     8 Params"
        print "1       1 X-t        1 Model     5 Params"
        print "2       1 X-t        2 Models    8 Params"
        print "3       1 X-t        2 Models    5 Params"
        print "4       2 X-t        1 Model    16 Params"
        print "5       2 X-t        1 Model    13 Params"
        print "6       2 X-t        1 Model    10 Params"


    #TODO: Optimize me please
    def PCS(self, parsedObj, type_of_fit, parsedObj2="None"):
        x    = parsedObj.getAllXarray()
        y    = parsedObj.getAllYarray()
        z    = parsedObj.getAllZarray()
        meas = parsedObj.getAllMeasarray()
        met  = parsedObj.getMetalLoc()
        ax   = parsedObj.getAxialVVU()
        rh   = parsedObj.getRhombicVVU()
        a    = parsedObj.getAlpha()
        b    = parsedObj.getBeta()
        g    = parsedObj.getGamma()
        #NOTE: 1 X-t, 1 Model, 8 Params
        if type_of_fit == 0:
            #NOTE: Agrees with Numbat
            p0 = zeros(8)
            p0[0], p0[1], p0[2], p0[3], p0[4], p0[5], p0[6], p0[7] = \
                                         met[0], met[1], met[2], ax, rh, a, b ,g
            soln,cov,info,mesg,success = leastsq(PCS1M1S, p0, \
                                              args=(meas, x,y,z), full_output=1)
            qual = FitSummary(soln,cov,info,mesg,success,p0, meas, type_of_fit)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln)
        #NOTE: 1 X-t, 1 Model, 5 Params (metal position fixed)
        if type_of_fit == 1:
            p0 = zeros(5)
            p0[0], p0[1], p0[2], p0[3], p0[4]  = ax, rh, a, b ,g
            soln,cov,info,mesg,success = leastsq(PCS1M1SFM, p0, \
                      args=(meas, met[0], met[1], met[2], x,y,z), full_output=1)
            qual = FitSummary(soln,cov,info,mesg,success,p0, meas, type_of_fit)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln)
        #NOTE: 1 X-t, 2 Models  , 8 Params (DIMER fit)
        if type_of_fit == 2:
            #NOTE: Agrees with Numbat
            x2 = parsedObj2.getAllXarray()
            y2 = parsedObj2.getAllYarray()
            z2 = parsedObj2.getAllZarray()
            p0 = zeros(8)
            p0[0], p0[1], p0[2], p0[3], p0[4], p0[5], p0[6], p0[7] = \
                                         met[0], met[1], met[2], ax, rh, a, b ,g
            soln,cov,info,mesg,success = leastsq(PCS2M1S, p0, \
                                    args=(meas, x,y,z, x2,y2,z2), full_output=1)
            qual = FitSummary(soln,cov,info,mesg,success,p0, meas, type_of_fit)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parsedObj2)
        #NOTE: 1 X-t, 2 Models, 5 Params (metal position fixed, DIMER fit)
        if type_of_fit == 3:
            #NOTE: Agrees with Numbat
            x2 = parsedObj2.getAllXarray()
            y2 = parsedObj2.getAllYarray()
            z2 = parsedObj2.getAllZarray()
            p0 = zeros(5)
            p0[0], p0[1], p0[2], p0[3], p0[4] = ax, rh, a, b ,g
            soln,cov,info,mesg,success = leastsq(PCS2M1SFM, p0, args= \
                 (meas, met[0], met[1], met[2], x,y,z, x2,y2,z2), full_output=1)
            qual = FitSummary(soln,cov,info,mesg,success,p0, meas, type_of_fit)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parsedObj2)
        #NOTE: 2 X-t, 1 Model, 16 Params
        if type_of_fit == 4:
            #NOTE: Agrees with Numbat2Tens
            #NOTE: Adjusted max number of function evaluations
            met2 = parsedObj2.getMetalLoc()
            ax2  = parsedObj2.getAxialVVU()
            rh2  = parsedObj2.getRhombicVVU()
            a2   = parsedObj2.getAlpha()
            b2   = parsedObj2.getBeta()
            g2   = parsedObj2.getGamma()
            p0 = zeros(16)
            p0[0], p0[1],  p0[2],  p0[3],  p0[4],  p0[5],  p0[6], p0[7]  = \
             met[0], met[1], met[2], ax, rh, a, b ,g
            p0[8], p0[9], p0[10], p0[11], p0[12], p0[13], p0[14], p0[15] = \
             met2[0], met2[1], met2[2], ax2, rh2, a2, b2 ,g2
            soln,cov,info,mesg,success = leastsq(PCS1M2S, p0, args= \
                                    (meas, x,y,z), full_output=1, maxfev=100000)
            qual = FitSummary(soln,cov,info,mesg,success,p0, meas, type_of_fit)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parsedObj2)
        #NOTE: 2 X-t, 1 Model, 13 Params (metal site 1 is fixed)
        if type_of_fit == 5:
            #TODO: Compare results to Numbat2Tens
            #NOTE: Adjusted max number of function evaluations
            met2 = parsedObj2.getMetalLoc()
            ax2  = parsedObj2.getAxialVVU()
            rh2  = parsedObj2.getRhombicVVU()
            a2   = parsedObj2.getAlpha()
            b2   = parsedObj2.getBeta()
            g2   = parsedObj2.getGamma()
            p0 = zeros(13)
            p0[0], p0[1], p0[2], p0[3], p0[4], = ax, rh, a, b ,g
            p0[5], p0[6], p0[7], p0[8], p0[9], p0[10], p0[11], p0[12] = \
             met2[0], met2[1], met2[2], ax2, rh2, a2, b2 ,g2
            soln,cov,info,mesg,success = leastsq(PCS1M2SFMX1, p0, args= \
            (meas, met[0], met[1], met[2], x,y,z), full_output=1, maxfev=100000)
            qual = FitSummary(soln,cov,info,mesg,success,p0, meas, type_of_fit)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parsedObj2)
        #NOTE: 2 X-t, 1 Model, 10 Params (metal sites 1 and 2 are fixed)
        if type_of_fit == 6:
            met2 = parsedObj2.getMetalLoc()
            ax2  = parsedObj2.getAxialVVU()
            rh2  = parsedObj2.getRhombicVVU()
            a2   = parsedObj2.getAlpha()
            b2   = parsedObj2.getBeta()
            g2   = parsedObj2.getGamma()
            p0 = zeros(10)
            p0[0], p0[1], p0[2], p0[3], p0[4] = ax , rh,  a , b  ,g
            p0[5], p0[6], p0[7], p0[8], p0[9] = ax2, rh2, a2, b2 ,g2
            soln,cov,info,mesg,success = leastsq(PCS1M2SFMX2, p0, args= \
             (meas, met[0], met[1], met[2], met2[0], met2[1], met2[2], x,y,z), \
              full_output=1, maxfev=100000)
            qual = FitSummary(soln,cov,info,mesg,success,p0, meas, type_of_fit)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parsedObj2)

        if success == 1:
            return soln, qual
        else:
            print
            print "There were problems with the optimization."
            print "This could mean optimization failed."
            print "If in dobt, re-optimize with different initial conditions!"
            print
            return soln, qual



    def PREOptions(self):
        """
        Prints information on the types of fits that can be achieved from PCS
        data
        """
        print "0       1 centre     1 Model     Fixed constant      3 Params"
        print "1       1 centre     1 Model     Optimize constant   4 Params"
        print "2       1 centre     2 Models    Fixed constant      3 Params"
        print "3       1 X-t        2 Models    Optimize constant   4 Params"
        print "4       2 centres    1 Model     Fixed constant      6 Params"
        print "5       2 centres    1 Model     Optimize constant   7 Params"


    #TODO: Optimize me please
    def PRE(self, parsedObj, type_of_fit, parsedObj2="None"):
        x    = parsedObj.getAllXarray()
        y    = parsedObj.getAllYarray()
        z    = parsedObj.getAllZarray()
        meas = parsedObj.getAllMeasarray()
        c    = parsedObj.getConstant()
        #NOTE: Change below
        tol = parsedObj.getAllTolarray()
        #NOTE: 1 centre, 1 Model, Fixed constant, 3 Params
        if type_of_fit == 0:
            #NOTE: Agrees with PREfit (by fixing c to free optimized value.
            p0   = parsedObj.getMetalLoc()
            soln,cov,info,mesg,success = leastsq(PRE1M1SFC, p0, args= \
                                                 (meas, c,x,y,z, tol), full_output=1)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln)
        #NOTE:1 centre, 1 Model, Optimize constant, 4 Params
        if type_of_fit == 1:
            #NOTE: Agrees with PREfit
            met   = parsedObj.getMetalLoc()
            p0 = zeros(4)
            p0[0], p0[1], p0[2], p0[3] = met[0],met[1],met[2], c[0]
            soln,cov,info,mesg,success = leastsq(PRE1M1SOC, p0, args= \
                                                   (meas, x,y,z), full_output=1)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln)

        if type_of_fit == 2:
            #NOTE: Agrees with PREfit
            x2    = parsedObj2.getAllXarray()
            y2    = parsedObj2.getAllYarray()
            z2    = parsedObj2.getAllZarray()
            p0    = parsedObj2.getMetalLoc()
            soln,cov,info,mesg,success = leastsq(PRE2M1SFC, p0, args= \
                                        (meas, c,x,y,z,x2,y2,z2, tol), full_output=1)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parsedObj2)
        if type_of_fit == 3:
            #NOTE: Agrees with PREfit
            x2    = parsedObj2.getAllXarray()
            y2    = parsedObj2.getAllYarray()
            z2    = parsedObj2.getAllZarray()
            met   = parsedObj.getMetalLoc()
            p0 = zeros(4)
            p0[0], p0[1], p0[2], p0[3] = met[0],met[1],met[2], c[0]
            soln,cov,info,mesg,success = leastsq(PRE2M1SOC, p0, args= \
                                          (meas, x,y,z,x2,y2,z2), full_output=1)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parsedObj2)
        if type_of_fit == 4:
            p0 = zeros(6)
            m1 = parsedObj.getMetalLoc()
            m2 = parsedObj2.getMetalLoc()
            p0[0], p0[1], p0[2], p0[3], p0[4], p0[5]  = m1[0],m1[1],m1[2], \
                                                        m2[0],m2[1],m2[2]
            soln,cov,info,mesg,success = leastsq(PRE1M2SFC, p0, args= \
                                                 (meas, c,x,y,z, tol), full_output=1)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parsedObj2)
        if type_of_fit == 5:
            #NOTE: Agrees with PREFit.
            #NOTE: Optimization may not converge to same params (PREfit/pyPara)
            p0 = zeros(7)
            m1 = parsedObj.getMetalLoc()
            m2 = parsedObj2.getMetalLoc()
            p0[0],p0[1],p0[2],p0[3],p0[4],p0[5],p0[6]  = m1[0],m1[1],m1[2], \
                                                   m2[0],m2[1],m2[2], c[0]
            soln,cov,info,mesg,success = leastsq(PRE1M2SOC, p0, args= \
            (meas, x,y,z), full_output=1)
            self.__UpdateFromFit(parsedObj, type_of_fit, soln, parseObj2)

#        if success == 1:
        return soln,cov,info,mesg,success
#        else:
#            msg1 = "Optimization failed. Please try again with different "
#            msg2 = "different intial starting guesses"
#            print msg1+msg2


    def RDCOptions(self):
        """
        Prints information on the types of fits that can be achieved from RDC
        data
        """
        print "0       1 A-t        1 Model     5 Params"


    #TODO: Optimize me please
    def RDC(self, parsedObj, type_of_fit, scal, parsedObj2="None"):
        x1,x2 = parsedObj.getAllXarray()
        y1,y2 = parsedObj.getAllYarray()
        z1,z2 = parsedObj.getAllZarray()
        meas  = parsedObj.getAllMeasarray()
        Dax   = parsedObj.getAxial()
        Drh   = parsedObj.getRhombic()
        a     = parsedObj.getAlpha()
        b     = parsedObj.getBeta()
        g     = parsedObj.getGamma()
        #NOTE: 1 A-t, 1 Model, 5 Params
        if type_of_fit == 0:
            p0 = zeros(5)
            p0[0], p0[1], p0[2], p0[3], p0[4] = Dax, Drh, a, b ,g
            soln,cov,info,mesg,success = leastsq(RDC1M1S, p0, \
                                              args=(meas, x1,y1,z1, \
                                                          x2,y2,z2, scal),\
                                                          full_output=1)
        print soln



    def PCS_PRE(self, parsedObj, type_of_fit, parsedObj2, parsedObj3, \
                    wt_pcs, wt_pre):
        x1       = parsedObj.getAllXarray()
        y1       = parsedObj.getAllYarray()
        z1       = parsedObj.getAllZarray()
        meas_pcs = parsedObj.getAllMeasarray()
        met1     = parsedObj.getMetalLoc()
        ax       = parsedObj.getAxialVVU()
        rh       = parsedObj.getRhombicVVU()
        a        = parsedObj.getAlpha()
        b        = parsedObj.getBeta()
        g        = parsedObj.getGamma()
        #Fit two sites to a a single model using PCS and PRE data Fixed C
        if type_of_fit == 0:
            x2       = parsedObj2.getAllXarray()
            y2       = parsedObj2.getAllYarray()
            z2       = parsedObj2.getAllZarray()
            meas_pre = parsedObj2.getAllMeasarray()
            c_val    = parsedObj2.getConstant()
            met2     = parsedObj3.getMetalLoc()
            ax2      = parsedObj3.getAxialVVU()
            rh2      = parsedObj3.getRhombicVVU()
            a2       = parsedObj3.getAlpha()
            b2       = parsedObj3.getBeta()
            g2       = parsedObj3.getGamma()

            p0 = zeros(16)
            p0[0], p0[1],  p0[2],  p0[3],  p0[4],  p0[5],  p0[6], p0[7]  = \
             met1[0], met1[1], met1[2], ax, rh, a, b ,g
            p0[8], p0[9], p0[10], p0[11], p0[12], p0[13], p0[14], p0[15] = \
             met2[0], met2[1], met2[2], ax2, rh2, a2, b2 ,g2

            soln,cov,info,mesg,success = leastsq(PCSPRE1M2SFC, p0, args= \
              (meas_pcs, meas_pre, x1,y1,z1, x2,y2,z2, c_val, wt_pcs, wt_pre), \
               full_output=1, maxfev=10000)
            #Note lowered fcalls
            #qual = FitSummary(soln,cov,info,mesg,success,p0, meas, type_of_fit)
            #self.__UpdateFromFit(parsedObj, soln, parsedObj2)
            return soln



    def PCS_RDC(self, parsedObj, type_of_fit):
        # NOTE: This method will come in future (when RDC calc is known to work
        pass


    def PRE_RDC(self, parsedObj, type_of_fit):
        # NOTE: This method will come in future (when RDC calc is known to work
        pass


    def __UpdateFromFit(self, parsedObj, type_of_fit,fitRes, parsedObj2="None"):
        """
        Update a paramagnetic object with post optimization values
        @param parsedObj  : A paramagnetic object
        @param type_of_fit: Type of fit performed previously
        @param fitRes     : The solution array post optimization
        @param parsedObj2 : [OPTIONAL] A second parsedObj
        """
        if parsedObj.getDataType() == 'pcs':
            if type_of_fit  == 0:
                ms = zeros(3)
                ms[0], ms[1], ms[2] = fitRes[0], fitRes[1], fitRes[2]
                parsedObj.setMetalLoc(ms)
                parsedObj.setAxial(FromVVU(fitRes[3]))
                parsedObj.setRhombic(FromVVU(fitRes[4]))
                parsedObj.setAlpha(FixAngle(fitRes[5]))
                parsedObj.setBeta(FixAngle(fitRes[6]))
                parsedObj.setGamma(FixAngle(fitRes[7]))
                self.__AnglesUTR(parsedObj)
            if type_of_fit == 1:
                parsedObj.setAxial(FromVVU(fitRes[0]))
                parsedObj.setRhombic(FromVVU(fitRes[1]))
                parsedObj.setAlpha(FixAngle(fitRes[2]))
                parsedObj.setBeta(FixAngle(fitRes[3]))
                parsedObj.setGamma(FixAngle(fitRes[4]))
                self.__AnglesUTR(parsedObj)
            if type_of_fit == 2:
                ms = zeros(3)
                ms[0], ms[1], ms[2] = fitRes[0], fitRes[1], fitRes[2]
                parsedObj.setMetalLoc(ms)
                parsedObj.setAxial(FromVVU(fitRes[3]))
                parsedObj.setRhombic(FromVVU(fitRes[4]))
                parsedObj.setAlpha(FixAngle(fitRes[5]))
                parsedObj.setBeta(FixAngle(fitRes[6]))
                parsedObj.setGamma(FixAngle(fitRes[7]))
                parsedObj2.setMetalLoc(ms)
                parsedObj2.setAxial(FromVVU(fitRes[3]))
                parsedObj2.setRhombic(FromVVU(fitRes[4]))
                parsedObj2.setAlpha(FixAngle(fitRes[5]))
                parsedObj2.setBeta(FixAngle(fitRes[6]))
                parsedObj2.setGamma(FixAngle(fitRes[7]))
                self.__AnglesUTR(parsedObj)
                self.__AnglesUTR(parsedObj2)
            if type_of_fit == 3:
                parsedObj.setAxial(FromVVU(fitRes[0]))
                parsedObj.setRhombic(FromVVU(fitRes[1]))
                parsedObj.setAlpha(FixAngle(fitRes[2]))
                parsedObj.setBeta(FixAngle(fitRes[3]))
                parsedObj.setGamma(FixAngle(fitRes[4]))
                parsedObj2.setAxial(FromVVU(fitRes[0]))
                parsedObj2.setRhombic(FromVVU(fitRes[1]))
                parsedObj2.setAlpha(FixAngle(fitRes[2]))
                parsedObj2.setBeta(FixAngle(fitRes[3]))
                parsedObj2.setGamma(FixAngle(fitRes[4]))
                self.__AnglesUTR(parsedObj)
                self.__AnglesUTR(parsedObj2)
            if type_of_fit == 4:
                ms1 = zeros(3)
                ms1[0], ms1[1], ms1[2] = fitRes[0], fitRes[1], fitRes[2]
                parsedObj.setMetalLoc(ms1)
                parsedObj.setAxial(FromVVU(fitRes[3]))
                parsedObj.setRhombic(FromVVU(fitRes[4]))
                parsedObj.setAlpha(FixAngle(fitRes[5]))
                parsedObj.setBeta(FixAngle(fitRes[6]))
                parsedObj.setGamma(FixAngle(fitRes[7]))
                ms2 = zeros(3)
                ms2[0], ms2[1], ms2[2] = fitRes[8], fitRes[9], fitRes[10]
                parsedObj2.setMetalLoc(ms2)
                parsedObj2.setAxial(FromVVU(fitRes[11]))
                parsedObj2.setRhombic(FromVVU(fitRes[12]))
                parsedObj2.setAlpha(FixAngle(fitRes[13]))
                parsedObj2.setBeta(FixAngle(fitRes[14]))
                parsedObj2.setGamma(FixAngle(fitRes[15]))
                self.__AnglesUTR(parsedObj)
                self.__AnglesUTR(parsedObj2)

        if parsedObj.getDataType() == 'pre':
            if type_of_fit == 0:
                ms = zeros(3)
                ms[0], ms[1], ms[2] = fitRes[0], fitRes[1], fitRes[2]
                parsedObj.setMetalLoc(ms)
            if type_of_fit == 1:
                ms[0], ms[1], ms[2] = fitRes[0], fitRes[1], fitRes[2]
                parsedObj.setMetalLoc(ms)
                parsedObj.setConstant(fitRes[3])
            if type_of_fit == 2:
                ms = zeros(3)
                ms[0], ms[1], ms[2] = fitRes[0], fitRes[1], fitRes[2]
                parsedObj.setMetalLoc(ms)
                parsedObj2.setMetalLoc(ms)
            if type_of_fit == 3:
                ms = zeros(3)
                ms[0], ms[1], ms[2] = fitRes[0], fitRes[1], fitRes[2]
                parsedObj.setMetalLoc(ms)
                parsedObj.setConstant(fitRes[3])
                parsedObj2.setMetalLoc(ms)
                parsedObj2.setConstant(fitRes[3])
            if type_of_fit == 4:
                ms1 = zeros(3)
                ms2 = zeros(3)
                ms1[0], ms1[1], ms1[2] = fitRes[0], fitRes[1], fitRes[2]
                ms2[0], ms2[1], ms2[2] = fitRes[3], fitRes[4], fitRes[5]
                parsedObj.setMetalLoc(ms1)
                parsedObj2.setMetalLoc(ms2)
            if type_of_fit == 5:
                ms1 = zeros(3)
                ms2 = zeros(3)
                ms1[0], ms1[1], ms1[2] = fitRes[0], fitRes[1], fitRes[2]
                ms2[0], ms2[1], ms2[2] = fitRes[3], fitRes[4], fitRes[5]
                parsedObj.setMetalLoc(ms1)
                parsedObj.setConstant(fitRes[6])
                parsedObj2.setMetalLoc(ms2)
                parsedObj2.setConstant(fitRes[6])

        if parsedObj.getDataType() == 'rdc':
            #TODO: Finish/Add this
            pass



    def __AnglesUTR(self, parsedObj):
        a = parsedObj.getAlpha()
        b = parsedObj.getBeta()
        g = parsedObj.getGamma()
        Dx = -ToVVU(parsedObj.getAxial())/3.0 + \
              ToVVU(parsedObj.getRhombic())/2.0
        Dy = -ToVVU(parsedObj.getAxial())/3.0 - \
              ToVVU(parsedObj.getRhombic())/2.0
        Dz = 2.0/3.0*ToVVU(parsedObj.getAxial())
        aDx, aDy, aDz = abs(Dx), abs(Dy), abs(Dz)

        if (aDz >= aDy) and (aDy >= aDx):
            print "UTR Case1"
        if (aDz >= aDx)and (aDx >= aDy):
            g = g + 90.0
            Dy, Dx = SwapVals(Dy,Dx)
            print "UTR Case2"
        if (aDy >= aDz) and (aDz >= aDx):
            Dy, Dz = SwapVals(Dy,Dz)
            rX90 = RotX90()
            rZYZ = ZYZRot(a,b,g)
            nR = mat(rX90) * mat(rZYZ)
            a,b,g = ABGFromRotMatrixZYZ(nR)
            a,b,g = math.degrees(a), math.degrees(b), math.degrees(g)
            print "UTR Case3"
        if (aDy >= aDx) and (aDx >= aDz):
            g = g + 90.0
            Dy, Dx = SwapVals(Dy, Dx)
            Dz, Dx = SwapVals(Dz, Dx)
            rY90 = RotY90()
            rZYZ = ZYZRot(a,b,g)
            nR = mat(rY90) * mat(rZYZ)
            a,b,g = ABGFromRotMatrixZYZ(nR)
            a,b,g = math.degrees(a), math.degrees(b), math.degrees(g)
            print "UTR Case4"
        if(aDx >= aDz) and (aDz >= aDy):
            g = g + 90.0
            Dy, Dx = SwapVals(Dy, Dx)
            Dy, Dz = SwapVals(Dy, Dz)
            rX90 = RotX90()
            rZYZ = ZYZRot(a,b,g)
            nR = mat(rX90) * mat(rZYZ)
            a,b,g = ABGFromRotMatrixZYZ(nR)
            a,b,g = math.degrees(a), math.degrees(b), math.degrees(g)
            print "UTR Case5"
        if(aDx >= aDy) and (aDy >= aDz):
            Dz, Dx = SwapVals(Dz, Dx)
            rY90 = RotY90()
            rZYZ = ZYZRot(a,b,g)
            nR =  mat(rY90)* mat(rZYZ)
            a,b,g = ABGFromRotMatrixZYZ(nR)
            a,b,g = math.degrees(a), math.degrees(b), math.degrees(g)
            print "UTR Case6"

        #Axial and Rhombic are now in UTR
        Ax = Dz - (Dx + Dy)/2.0;
        Rh = Dx - Dy;
        Ax, Rh = FromVVU(Ax), FromVVU(Rh)

        # Make Euler angles in 0-360 after manipulation.
        a = FixAngle(a)
        b = FixAngle(b)
        g = FixAngle(g)

        # Do manipulations such that A,B,G in 0-180
        if a >= 0.0 and a < 180.0:                              # a in 0,   180
            if b >= 0.0  and  b < 180.0:                        # b in 0,   180
                if g >= 0.0  and  g < 180.0:                    # g in 0,   180
                    pass
                else:                                           # g in 180, 360
                    g = g + 180.0
            else:                                               # b in 180, 360
                if g >= 0.0 and g < 180.0:                      # g in 0,   180
                    b =  b + 180.0
                    g = -g +180
                else:                                           # g in 180, 360
                    b = b + 180.0
                    g = -g
        else:                                                   # a in 180, 360
            if b >= 0 and  b < 180.0:                           # b in 0,   180
                if g >=0  and  g < 180.0:                       # g in 0,   180
                    a =  a + 180.0
                    b = -b + 180.0
                    g = -g + 180.0
                else:                                           # g in 180, 360
                    a =  a + 180.0
                    b = -b + 180.0
                    g = -g
            else:                                               # b in 180, 360
                if g >= 0 and  g < 180.0:                       # g in 0,   180
                    a =  a + 180.0
                    b = -b
                    g =  g
                else:                                           # g in 180, 360
                    a = a + 180.0
                    b = -b
                    g = g + 180.0

        # Important. Fix to 0-360 to get in UTR (really 0-180).
        a = FixAngle(a)
        b = FixAngle(b)
        g = FixAngle(g)

        #Update for UTR!
        parsedObj.setAlpha(a)
        parsedObj.setBeta(b)
        parsedObj.setGamma(g)
        parsedObj.setAxial(Ax)
        parsedObj.setRhombic(Rh)


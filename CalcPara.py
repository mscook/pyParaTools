#TODO: Test and check current test coverage
#TODO: Check that the calculation over multiple ParsedObjs works correctly
#TODO: Ensure that the average calculated is true
#NOTE: All above can be adressed in the testing framework test_CalcPara.py

import math
import sys
from   ParaUtils import *


class CalcPara:
    """Calculates paramagnetic observables, PCS, PRE or RDC
       given a set of parameters"""


    def PCS(self, ParsedObj, convention, ParsedObjL=[]):
        """
        Calculate the PCS for a given set of spins/parameters in a given
        convention
        @param ParsedObj:  A container of parsed dataset/structure/params
        @type  ParsedObj:  ParaParser object
        @param convention: The X-tensor convention (ZXZ or ZYZ)
        @type  convention: string
        @param ParsedObjL: [OPTIONAL] A list of ParaParser objects. When
                           calculating the PCS on many ParaParser objects
        @type  ParsedObjL: list
        """
        ParsedObjL.insert(0,ParsedObj)
        average = []
        for i in range(0, len(ParsedObjL)):
            if convention == 'ZXZ':
                rot  = ZXZRot(ParsedObjL[i].getAlpha(), \
                              ParsedObjL[i].getBeta(),  \
                              ParsedObjL[i].getGamma())
            elif convention == 'ZYZ':
                rot  = ZYZRot(ParsedObjL[i].getAlpha(), \
                              ParsedObjL[i].getBeta(),  \
                              ParsedObjL[i].getGamma())
            else:
                print "Non-supported convention selected in PCS calculation"
                print "Exiting"
                sys.exit(0)
            metalPos = ParsedObjL[i].getMetalLoc()
            Xax      = ParsedObjL[i].getAxialVVU()
            Xrh      = ParsedObjL[i].getRhombicVVU()
            data     = ParsedObjL[i].getParsed()
            for spin in range(0, len(data)):
                cur = data[spin].getCoord()
                X   = cur[0] - metalPos[0]
                Y   = cur[1] - metalPos[1]
                Z   = cur[2] - metalPos[2]
                x_t = rot[0][0]*X + rot[0][1]*Y + rot[0][2]*Z
                y_t = rot[1][0]*X + rot[1][1]*Y + rot[1][2]*Z
                z_t = rot[2][0]*X + rot[2][1]*Y + rot[2][2]*Z
                r2  = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
                r5  = (r2*r2) * math.sqrt(r2)
                tmp = 1.0/r5
                pcs = tmp*(Xax*(3.0*z_t*z_t-r2)+Xrh*1.5*(x_t*x_t-y_t*y_t))
                data[spin].setCVal(pcs)
                if i == 0:
                    average.append(pcs)
                else:
                    average[spin] = average[spin] + pcs
        return average


    def PRE(self, ParsedObj, ParsedObjL=[]):
        """
        Calculate the PRE for a given set of spins/parameters
        @param ParsedObj:  A container of parsed dataset/structure/params
        @type  ParsedObj:  ParaParser object
        @param ParsedObjL: [OPTIONAL] A list of ParsedObj
        @type  ParsedObjL: list
        """
        ParsedObjL.insert(0,ParsedObj)
        average = []
        for i in range(0, len(ParsedObjL)):
            metalPos = ParsedObjL[i].getMetalLoc()
            data     = ParsedObjL[i].getParsed()
            for spin in range(0, len(data)):
                cur = data[spin].getCoord()
                r2  = (metalPos[0] - cur[0])**2 + (metalPos[1] -
                       cur[1])**2 + (metalPos[2] - cur[2])**2
                r   = math.sqrt(r2)
                pre = ParsedObjL[i].getConstant()/(r**6)
                data[spin].setCVal(pre)
                if i == 0:
                    average.append(pre)
                else:
                    average[spin] = average[spin] + pre
        return average


    def RDC(self, ParsedObj, convention, exp_type='HN', Nmgr_type=1,
                ParsedObjL=[]):
        """
        Calculate the RDCs for a given set of spins/parameters in a given
        convention
        @param ParsedObj:  A container of parsed dataset/structure/params
        @type  ParsedObj:  ParaParser object
        @param convention: The Alignment-tensor convention (ZXZ or ZYZ)
        @type  convention:.string
        @param exp_type:   [DEFAULT] The type of RDC experiment: NH coupled RDC
        @type  exp_type:   string
        @param Nmgr_type:  [DEFAULT]The Nitrogen isotope: N15 isotope labelling
        @type  Nmgr_type:..integer
        @param ParsedObjL: [OPTIONAL] A list of ParsedObj
        @type  ParsedObjL: list
        """
        ParsedObjL.insert(0,ParsedObj)
        average = []
        for i in range(0, len(ParsedObjL)):
            if convention == 'ZXZ':
                rot  = ZXZRot(ParsedObjL[i].getAlpha(), ParsedObjL[i].getBeta(),
                    ParsedObjL[i].getGamma())
            elif convention == 'ZYZ':
                rot  = ZYZRot(ParsedObjL[i].getAlpha(), ParsedObjL[i].getBeta(),
                    ParsedObjL[i].getGamma())
            else:
                print "Non-supported convention selected in RDC calculation"
                print "Exiting"
                sys.exit(0)
            if exp_type != 'HN':
                print "Non-supported RDC experiment type"
                print "Exiting"
                sys.exit(0)
            Dax   = ParsedObjL[i].getAxial()
            Drh   = ParsedObjL[i].getRhombic()
            data  = ParsedObjL[i].getParsed()
            B0    = ParsedObjL[i].getB0()
            temp  = ParsedObjL[i].getTemp()
            S     = ParsedObjL[i].getOrder()
            g1    = lookupMGR(exp_type[0])
            g2    = lookupMGR(exp_type[1])[Nmgr_type]
            scal  = rdcScal(S, g1, g2, B0, temp)
            for spin in range(0, len(data)):
                cur1, cur2 = data[spin].getCoord()
                X    = cur1[0] - cur2[0]
                Y    = cur1[1] - cur2[1]
                Z    = cur1[2] - cur2[2]
                x_t  = rot[0][0]*X + rot[0][1]*Y +rot[0][2]*Z
                y_t  = rot[1][0]*X + rot[1][1]*Y +rot[1][2]*Z
                z_t  = rot[2][0]*X + rot[2][1]*Y +rot[2][2]*Z
                r2   = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
                r5   = (r2*r2) * math.sqrt(r2)
                tmp  = 1.0/r5
                urdc = tmp*(Dax*(3.0*z_t*z_t-r2)+Drh*1.5*(x_t*x_t-y_t*y_t))
                rdc  = urdc*scal
                data[spin].setCVal(rdc)
                if i == 0:
                    average.append(rdc)
                else:
                    average[spin] = average[spin] + rdc
        return average


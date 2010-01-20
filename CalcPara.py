import math
import sys
from   ParaUtils import *


"""Calculates paramagnetic observables given a set of parameters"""
#TODO: Investigate using protected data ie _var - *DONE*
#TODO: Test and check current test coverage

class CalcPara:


    def PCS(self, parsed_pcs_data, convention):
        """
        Calculate the PCS for a given set of spins/parameters in a given
        convention
        @param parsed_pcs_data: A container of parsed dataset/structure/params
        @type parsed_pcs_data : ParaParser object
        @param convention:    : The X-tensor convention (ZXZ or ZYZ)
        @type convention:     : string
        """
        if convention == 'ZXZ':
            rot  = ZXZRot(parsed_pcs_data.getAlpha(), parsed_pcs_data.getBeta(),
                parsed_pcs_data.getGamma())
        elif convention == 'ZYZ':
            rot  = ZYZRot(parsed_pcs_data.getAlpha(), parsed_pcs_data.getBeta(),
                parsed_pcs_data.getGamma())
        else:
            print "Non-supported convention selected in PCS calculation"
            print "Exiting"
            sys.exit(0)
        metalPos = parsed_pcs_data.getMetalLoc()
        Xax    = parsed_pcs_data.getAxialVVU()
        Xrh  = parsed_pcs_data.getRhombicVVU()
        data     = parsed_pcs_data.getParsed()
        for pObject in range(0, len(data)):
            cur = data[pObject].getCoord()
            X   = cur[0] - metalPos[0]
            Y   = cur[1] - metalPos[1]
            Z   = cur[2] - metalPos[2]
            x_t = rot[0][0]*X + rot[0][1]*Y + rot[0][2]*Z
            y_t = rot[1][0]*X + rot[1][1]*Y + rot[1][2]*Z
            z_t = rot[2][0]*X + rot[2][1]*Y + rot[2][2]*Z
            r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
            r5 = (r2*r2) * math.sqrt(r2)
            tmp = 1.0/r5
            pcs = tmp*(Xax * (3.0*z_t*z_t -r2) + Xrh*1.5*(x_t*x_t - y_t*y_t))
            data[pObject].setCVal(pcs)


    def PRE(self, parsed_pre_data):
        """
        Calculate the PREfor a given set of spins/parameters
        @param parsed_pre_data: A container of parsed dataset/structure/params
        @type parsed_pre_data:  ParaParser object
        """
        metalPos = parsed_pre_data.getMetalLoc()
        data     = parsed_pre_data.getParsed()
        for pObject in range(0, len(data)):
            cur = data[pObject].getCoord()
            r2 = (metalPos[0] - cur[0])**2 + (metalPos[1] -
                cur[1])**2 + (metalPos[2] - cur[2])**2
            r = math.sqrt(r2)
            pre = parsed_pre_data.getConstant()/(r**6)
            data[pObject].setCVal(pre)


    #FIXME: The RDC method is not finished/working. Must check!
    def RDC(parsed_rdc_data, convention, B0, temp, Nmgr_type=1):
        """
        FIXME
        @param parsed_rdc_data: A container of parsed dataset/structure/params
        @type parsed_rdc_data:  ParaParser object
        @param convention:      The X-tensor convention (ZXZ or ZYZ)
        @type convention:.......string
        @param B0:..............The magnetic field strength (Tesla)
        @type B0:               float
        @param temp:            The temperature (Kelvin)
        @type temp:             float
        @param Nmgr_type=1:     The Nitrogen isotope (N14 or N15)
        @type Nmgr_type=1:......boolean
        """
        if convention == 'ZXZ':
            rot  = ZXZRot(parsed_rdc_data.getAlpha(), parsed_rdc_data.getBeta(),
                parsed_rdc_data.getGamma())
        elif convention == 'ZYZ':
            rot  = ZYZRot(parsed_rdc_data.getAlpha(), parsed_rdc_data.getBeta(),
                parsed_rdc_data.getGamma())
        else:
            print "Non-supported convention selected in RDC calculation"
            print "Exiting"
            sys.exit(0)
        Dax      = parsed_rdc_data.getAxial()
        Drh      = parsed_rdc_data.getRhombic()
        data     = parsed_rdc_data.getParsed()
        #FIXME: Add this code. Below. Currently set to 0.0
        gH       = 0.0
        gN       = 0.0
        for pObject in range(0, len(data)):
            cur1, cur2 = data[pObject].getCoord()
            X = cur1[0] - cur2[0]
            Y = cur1[1] - cur2[1]
            Z = cur1[2] - cur2[2]
            x_t = rot[0][0]*X + rot[0][1]*Y +rot[0][2]*Z
            y_t = rot[1][0]*X + rot[1][1]*Y +rot[1][2]*Z
            z_t = rot[2][0]*X + rot[2][1]*Y +rot[2][2]*Z
            r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
            r5 = (r2*r2) * math.sqrt(r2)
            tmp = 1.0/r2
            us_rdc= (Dax * (3.0*z_t*z_t -r2) + Drh*1.5*(x_t*x_t - y_t*y_t))*tmp
            scal = RDCScal(B0, temp, gH, gN)
            rdc = scal*us_rdc
            data[pObject].setCVal(rdc)


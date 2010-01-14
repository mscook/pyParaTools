import math

class CalcPara:

#    plankc = 1.05457148e-34     # is hbar
#    kboltz = 1.3806503e-23
#    distNH = 1.03e-10
#    vvu_s  = 3.76991118431e-35  # not exactly sure where this comes from...


    def ZXZRot(self, A,B,G):
        rot = [None]*3
        for i in range(3):
            rot[i] = [None] * 3
        ca = math.cos(math.radians(A))
        cb = math.cos(math.radians(B))
        cg = math.cos(math.radians(G))
        sa = math.sin(math.radians(A))
        sb = math.sin(math.radians(B))
        sg = math.sin(math.radians(G))
        rot[0][0] = ( cg * ca) - (cb * sa * sg)
        rot[0][1] = ( cg * sa) + (cb * ca * sg)
        rot[0][2] = ( sg * sb)
        rot[1][0] = (-sg * ca) - (cb * sa * cg)
        rot[1][1] = (-sg * sa) + (cb * ca * cg)
        rot[1][2] = ( cg * sb)
        rot[2][0] = ( sb * sa)
        rot[2][1] = (-sb * ca)
        rot[2][2] =   cb
        return rot


    def ZYZRot(self,A,B,G):
        rot = [None]*3
        for i in range(3):
            rot[i] = [None] * 3
        ca = math.cos(math.radians(A))
        cb = math.cos(math.radians(B))
        cg = math.cos(math.radians(G))
        sa = math.sin(math.radians(A))
        sb = math.sin(math.radians(B))
        sg = math.sin(math.radians(G))
        rot[0][0] = (-sg * sa) + (cb * ca * cg)
        rot[0][1] = ( sg * ca) + (cb * sa * cg)
        rot[0][2] = ( -cg * sb)
        rot[1][0] = (-cg * sa) - (cb * ca * sg)
        rot[1][1] = (cg * ca) - (cb * sa * sg)
        rot[1][2] = ( sg * sb)
        rot[2][0] = ( sb * ca)
        rot[2][1] = (sb * sa)
        rot[2][2] =   cb
        return rot


#    def RDCScal(self, B0, temp, gH, gN):
#        a = -1*(B0**2)*gammaH*gammaN*(plankc)
#        b = ((distNH**3)*kboltz*temp)*(120*math.pi**2)
#        return (a/b)*vvu_s


    def PCSZXZ(self, parsed_pcs_data):
        zxz_rot  = self.ZXZRot(parsed_pcs_data.getAlpha(), parsed_pcs_data.getBeta(), parsed_pcs_data.getGamma())
        metalPos = parsed_pcs_data.getMetalLoc()
        Xax    = parsed_pcs_data.getAxial()
        Xrh  = parsed_pcs_data.getRhombic()
        data     = parsed_pcs_data.getParsed()

        for pObject in range(0, len(data)):
            cur = data[pObject].getCoord()
            X   = cur[0] - metalPos[0]
            Y   = cur[1] - metalPos[1]
            Z   = cur[2] - metalPos[2]
            x_t = zxz_rot[0][0]*X + zxz_rot[0][1]*Y + zxz_rot[0][2]*Z
            y_t = zxz_rot[1][0]*X + zxz_rot[1][1]*Y + zxz_rot[1][2]*Z
            z_t = zxz_rot[2][0]*X + zxz_rot[2][1]*Y + zxz_rot[2][2]*Z
            r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
            r5 = (r2*r2) * math.sqrt(r2)
            tmp = 1.0/r5
            pcs = tmp*(Xax * (3.0*z_t*z_t -r2) + Xrh*1.5*(x_t*x_t - y_t*y_t))
            print pcs

    def PCSZYZ(self, parsed_pcs_data):
        zyz_rot  = self.ZYZRot(parsed_pcs_data.getAlpha(), parsed_pcs_data.getBeta(), parsed_pcs_data.getGamma())
        metalPos = parsed_pcs_data.getMetalLoc()
        Xax      = parsed_pcs_data.getAxial()
        Xrh      = parsed_pcs_data.getRhombic()
        data     = parsed_pcs_data.getParsed()

        for pObject in range(0, len(data)):
            cur = data[pObject].getCoord()
            X   = round(cur[0], 3) - metalPos[0]
            Y   = round(cur[1], 3) - metalPos[1]
            Z   = round(cur[2], 3) - metalPos[2]
            x_t = zyz_rot[0][0]*X + zyz_rot[0][1]*Y + zyz_rot[0][2]*Z
            y_t = zyz_rot[1][0]*X + zyz_rot[1][1]*Y + zyz_rot[1][2]*Z
            z_t = zyz_rot[2][0]*X + zyz_rot[2][1]*Y + zyz_rot[2][2]*Z
            r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
            r5 = (r2*r2) * math.sqrt(r2)
            tmp = 1.0/r5
            pcs = tmp*(Xax * (3.0*z_t*z_t -r2) + Xrh*1.5*(x_t*x_t - y_t*y_t))
            print data[pObject].getName(),data[pObject].getId(), pcs


#    def RDCZXZ(parsed_rdc_data, B0, temp, Nmgr_type=1):
#        zxz_rot  = self.ZXZRot(parsed_rdc_data.getAlpha(), parsed_rdc_data.getBeta(), parsed_rdc_data.getGamma())
#        Dax      = parsed_pcs_data.getAxial()
#        Drh      = parsed_pcs_data.getRhombic()
#        data     = parsed_pcs_data.getParsed()
#        gH       =
#        gN       =

#        for pObject in range(0, len(data)):
#            cur1, cur2 = data[pObject].getCoord()
#            X = cur1[0] - cur2[0]
#            Y = cur1[1] - cur2[1]
#            Z = cur1[2] - cur2[2]
#            x_t = zxz_rot[0][0]*X + zxz_rot[0][1]*Y +zxz_rot[0][2]*Z
#            y_t = zxz_rot[1][0]*X + zxz_rot[1][1]*Y +zxz_rot[1][2]*Z
#            z_t = zxz_rot[2][0]*X + zxz_rot[2][1]*Y +zxz_rot[2][2]*Z
#            r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
#            r5 = (r2*r2) * math.sqrt(r2)
#            tmp = 1.0/r2
#            us_rdc= (Dax * (3.0*z_t*z_t -r2) + Drh*1.5*(x_t*x_t - y_t*y_t))*tmp
#            scal = RDCScal(B0, temp, gH, gN):
#            rdc = scal*us_rdc
#            print rdc

#    def RDCZYZ(parsed_rdc_data):
#        zyz_rot  = self.ZYZRot(parsed_rdc_data.getAlpha(), parsed_rdc_data.getBeta(), parsed_rdc_data.getGamma())
#        Dax      = parsed_pcs_data.getAxial()
#        Drh      = parsed_pcs_data.getRhombic()
#        data     = parsed_pcs_data.getParsed()
#        gH       =
#        gN       =

#        for pObject in range(0, len(data)):
#            cur1, cur2 = data[pObject].getCoord()
#            X = cur1[0] - cur2[0]
#            Y = cur1[1] - cur2[1]
#            Z = cur1[2] - cur2[2]
#            x_t = zyz_rot[0][0]*X + zyz_rot[0][1]*Y +zyz_rot[0][2]*Z
#            y_t = zyz_rot[1][0]*X + zyz_rot[1][1]*Y +zyz_rot[1][2]*Z
#            z_t = zyz_rot[2][0]*X + zyz_rot[2][1]*Y +zyz_rot[2][2]*Z
#            r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
#            r5 = (r2*r2) * math.sqrt(r2)
#            tmp = 1.0/r2
#            us_rdc= (Dax * (3.0*z_t*z_t -r2) + Drh*1.5*(x_t*x_t - y_t*y_t))*tmp
#            scal = RDCScal(B0, temp, gH, gN):
#            rdc = scal*us_rdc
#            print rdc



    def PRE(self, parsed_pre_data):
        metalPos = parsed_pre_data.getMetalLoc()
        data     = parsed_pre_data.getParsed()
        for pObject in range(0, len(data)):
            cur = data[pObject].getCoord()
            r2 = (metalPos[0] - cur[0])**2 + (metalPos[1] - cur[1])**2 + (metalPos[2] - cur[2])**2
            r = math.sqrt(r2)
            print parsed_pre_data.getConstant()/(r**6)


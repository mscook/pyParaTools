import math

###############################################################################
#
# Calculate: 
#   * Pseudocontact Shift                   (PCS)
#   * Residual Dipolar Coupling             (RDC)
#   * Paramagnetic Relaxation Enhancement   (PRE)
# Authour: Mitchell Stanton-Cook
###############################################################################



###############################################################################
#
# Calculate the PCS given for ZXZ rotation
#
###############################################################################
def calcsPCSZXZ(angles, rH, rM, Dax, Drh, H):

	PI = math.pi

	rot = [["00", "01", "02"], ["10", "11", "12"], ["20", "21", "22"]] 
	
	ca = math.cos(angles[0]) ### FIXED !! 
	cb = math.cos(angles[1]) ### DON'T CONVERT TO RADIANS TWICE !!!
	cg = math.cos(angles[2])

	sa = math.sin(angles[0])
	sb = math.sin(angles[1])
	sg = math.sin(angles[2])

	rot[0][0] = ( cg * ca) - (cb * sa * sg)
	rot[0][1] = ( cg * sa) + (cb * ca * sg)
	rot[0][2] = ( sg * sb)
	rot[1][0] = (-sg * ca) - (cb * sa * cg)
	rot[1][1] = (-sg * sa) + (cb * ca * cg)
	rot[1][2] = ( cg * sb)
	rot[2][0] = ( sb * sa)
	rot[2][1] = (-sb * ca)
	rot[2][2] =   cb

	X = rH[0] - rM[0]
	Y = rH[1] - rM[1]
	Z = rH[2] - rM[2] 

	x_t = rot[0][0]*X + rot[0][1]*Y +rot[0][2]*Z
	y_t = rot[1][0]*X + rot[1][1]*Y +rot[1][2]*Z
	z_t = rot[2][0]*X + rot[2][1]*Y +rot[2][2]*Z

	r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
	r5 = (r2*r2) * math.sqrt(r2)

	tmp = 1.0/r5

	pcs = tmp*(Dax * (3.0*z_t*z_t -r2) + Drh*1.5*(x_t*x_t - y_t*y_t))
	return pcs



###############################################################################
#
# Calculate the PCS given for ZYZ rotation
#
###############################################################################
def calcsPCSZYZ(angles, rH, rM, Dax, Drh):

        PI = math.pi
	scal_const = (1./((12*PI))*10000)

	Dax = Dax*scal_const
	Drh = Drh*scal_const	

        rot = [["00", "01", "02"], ["10", "11", "12"], ["20", "21", "22"]]

        ca = math.cos(angles[0])
        cb = math.cos(angles[1])
        cg = math.cos(angles[2])

        sa = math.sin(angles[0])
        sb = math.sin(angles[1])
        sg = math.sin(angles[2])

        #	-sin(g)*sin(a)+cos(b)*cos(a)*cos(g)
	rot[0][0] = (-sg * sa) + (cb * ca * cg)
	#	sin(g)*cos(a)+cos(b)*sin(a)cos(g)
        rot[0][1] = ( sg * ca) + (cb * sa * cg)
        #	-cos(g)*sin(b)
	rot[0][2] = ( -cg * sb)
	#	-cos(g)*sin(a)-cos(b)*cos(a)*sin(g)
        rot[1][0] = (-cg * sa) - (cb * ca * sg)
	#	cos(g)*cos(a)-cos(b)*sin(a)*sin(g)
        rot[1][1] = (cg * ca) - (cb * sa * sg)
	#	sin(g)*sin(b)
        rot[1][2] = ( sg * sb)
	#	sin(b)*cos(a)
        rot[2][0] = ( sb * ca)
	#	sin(b)*sin(a)
        rot[2][1] = (sb * sa)
	#	cos(b)
        rot[2][2] =   cb

	X = rH[0] - rM[0]
        Y = rH[1] - rM[1]
        Z = rH[2] - rM[2]

        x_t = rot[0][0]*X + rot[0][1]*Y +rot[0][2]*Z
        y_t = rot[1][0]*X + rot[1][1]*Y +rot[1][2]*Z
        z_t = rot[2][0]*X + rot[2][1]*Y +rot[2][2]*Z

        r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
        r5 = (r2*r2) * math.sqrt(r2)

        tmp = 1.0/r5

        pcs = tmp*(Dax * (3.0*z_t*z_t -r2) + Drh*1.5*(x_t*x_t - y_t*y_t))
        return pcs



###############################################################################
#
# Calculate the PCS given the COMBINED rotation MATRIX
#
###############################################################################
def calcsPCS_MATRIX(rot_matrix, rH, rM, Dax, Drh, H):

	rot = [["00", "01", "02"], ["10", "11", "12"], ["20", "21", "22"]] 
	
	rot[0][0] = rot_matrix[0][0]
	rot[0][1] = rot_matrix[0][1] 
	rot[0][2] = rot_matrix[0][2] 
	rot[1][0] = rot_matrix[1][0]
	rot[1][1] = rot_matrix[1][1]
	rot[1][2] = rot_matrix[1][2]
	rot[2][0] = rot_matrix[2][0]
	rot[2][1] = rot_matrix[2][1]
	rot[2][2] = rot_matrix[2][2]


	X = rH[0] - rM[0]
	Y = rH[1] - rM[1]
	Z = rH[2] - rM[2] 

	x_t = rot[0][0]*X + rot[0][1]*Y +rot[0][2]*Z
	y_t = rot[1][0]*X + rot[1][1]*Y +rot[1][2]*Z
	z_t = rot[2][0]*X + rot[2][1]*Y +rot[2][2]*Z

	r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
	r5 = (r2*r2) * math.sqrt(r2)

	tmp = 1.0/r5

	pcs = tmp*(Dax * (3.0*z_t*z_t -r2) + Drh*1.5*(x_t*x_t - y_t*y_t))
	return pcs



###############################################################################
#
# Calculate the RDCs
#
###############################################################################
def calcRDCS(angles, rH, rN, Dax, Drh, H):

 #B0: 18.8
 #Magnetogyric ratio (rad/T s) H: 267.522e6
 #Magnetogyric ratio (rad/T s) N: 19.338e6
 #hbar (Plank constant/2Pi): 1.05457148e-34
 #Proton-Nitrogen distance:  1.03*e-10 [distHN]
 #Boltzman Constant: 1.3806503e-23 [kboltz]
 #T: 298
 #VVU 3.76991118431e-35
 
	PI = math.pi
	machB0 = 18.8
	temp =  298	
	gammaH = 267.522e6
	gammaN = 19.338e6
	plankc = 1.05457148e-34  # NOTE: HBAR
	kboltz = 1.3806503e-23
	distNH = 1.03e-10
        vvu_s = 3.76991118431e-35

	
	rot = [["00", "01", "02"], ["10", "11", "12"], ["20", "21", "22"]] 
	
	ca = math.cos(angles[0]) ### FIXED !! 
	cb = math.cos(angles[1]) ### DON'T CONVERT TO RADIANS TWICE !!!
	cg = math.cos(angles[2])

	sa = math.sin(angles[0])
	sb = math.sin(angles[1])
	sg = math.sin(angles[2])

	rot[0][0] = ( cg * ca) - (cb * sa * sg)
	rot[0][1] = ( cg * sa) + (cb * ca * sg)
	rot[0][2] = ( sg * sb)
	rot[1][0] = (-sg * ca) - (cb * sa * cg)
	rot[1][1] = (-sg * sa) + (cb * ca * cg)
	rot[1][2] = ( cg * sb)
	rot[2][0] = ( sb * sa)
	rot[2][1] = (-sb * ca)
	rot[2][2] =   cb

	X = rH[0] - rN[0]
	Y = rH[1] - rN[1]
	Z = rH[2] - rN[2] 

	x_t = rot[0][0]*X + rot[0][1]*Y +rot[0][2]*Z
	y_t = rot[1][0]*X + rot[1][1]*Y +rot[1][2]*Z
	z_t = rot[2][0]*X + rot[2][1]*Y +rot[2][2]*Z

	
	r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
	#r5 = (r2*r2) * math.sqrt(r2)
	tmp = 1.0/r2
	

	rdc_tmp = (Dax * (3.0*z_t*z_t -r2) + Drh*1.5*(x_t*x_t - y_t*y_t))*tmp
        a = -1*(machB0**2)*gammaH*gammaN*(plankc)
	b = ((distNH**3)*kboltz*temp)*(120*math.pi**2)
	rdc = (a/b)*vvu_s*rdc_tmp
	return rdc 

 
###############################################################################
#
# Calculate the RDCs given the COMBINED rotation MATRIX
#
###############################################################################
def calcRDCS_MATRIX(rot_matrix, rH, rN, Dax, Drh, H):

 #B0: 18.8
 #Magnetogyric ratio (rad/T s) H: 267.522e6
 #Magnetogyric ratio (rad/T s) N: 19.338e6
 #hbar (Plank constant/2Pi): 1.05457148e-34
 #Proton-Nitrogen distance:  1.03*e-10 [distHN]
 #Boltzman Constant: 1.3806503e-23 [kboltz]
 #T: 298
 #VVU 3.76991118431e-35
 
	PI = math.pi
	machB0 = 18.8
	temp =  298	
	gammaH = 267.522e6
	gammaN = 19.338e6
	plankc = 1.05457148e-34  # NOTE: HBAR
	kboltz = 1.3806503e-23
	distNH = 1.03e-10
        vvu_s = 3.76991118431e-35

	rot = [["00", "01", "02"], ["10", "11", "12"], ["20", "21", "22"]] 
	
	rot[0][0] = rot_matrix[0][0]
	rot[0][1] = rot_matrix[0][1] 
	rot[0][2] = rot_matrix[0][2] 
	rot[1][0] = rot_matrix[1][0]
	rot[1][1] = rot_matrix[1][1]
	rot[1][2] = rot_matrix[1][2]
	rot[2][0] = rot_matrix[2][0]
	rot[2][1] = rot_matrix[2][1]
	rot[2][2] = rot_matrix[2][2]


	X = rH[0] - rN[0]
	Y = rH[1] - rN[1]
	Z = rH[2] - rN[2] 

	x_t = rot[0][0]*X + rot[0][1]*Y +rot[0][2]*Z
	y_t = rot[1][0]*X + rot[1][1]*Y +rot[1][2]*Z
	z_t = rot[2][0]*X + rot[2][1]*Y +rot[2][2]*Z

	
	r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)
	#r5 = (r2*r2) * math.sqrt(r2)
	tmp = 1.0/r2
	

	rdc_tmp = (Dax * (3.0*z_t*z_t -r2) + Drh*1.5*(x_t*x_t - y_t*y_t))*tmp
        a = -1*(machB0**2)*gammaH*gammaN*(plankc)
	b = ((distNH**3)*kboltz*temp)*(120*math.pi**2)
	rdc = (a/b)*vvu_s*rdc_tmp
	return rdc 



###############################################################################
#
# Calculate the scaling of X-tensor to the Alignment Tensor
#
###############################################################################	
def RDC_scal(Dax, Drh):

 #B0: 18.8
 #Magnetogyric ratio (rad/T s) H: 267.522e6
 #Magnetogyric ratio (rad/T s) N: 19.338e6
 #hbar (Plank constant/2Pi): 1.05457148e-34
 #Proton-Nitrogen distance:  1.03*e-10 [distHN]
 #Boltzman Constant: 1.3806503e-23 [kboltz]
 #T: 298
 #VVU 3.76991118431e-35

        PI = math.pi
        machB0 = 18.8
        temp =  298
        gammaH = 267.522e6
        gammaN = 19.338e6
        plankc = 1.05457148e-34  # NOTE: HBAR
        kboltz = 1.3806503e-23
        distNH = 1.03e-10
        vvu_s = 3.76991118431e-35
	
	
        a = -1*(machB0**2)*gammaH*gammaN*(plankc)
        b = ((distNH**3)*kboltz*temp)*(120*math.pi**2)
        return (a/b)*vvu_s*Dax, (a/b)*vvu_s*Drh


def calc_PRE(c, rS, rM):
    """Calculate the PRE effect using c/r^6
        Inputs:
            * c the constant
            * rS, the x,y,z of the spin of interest
            * rM, the x,y,z location of the paramagnetic centre """
    r2 = (rM[0] - rS[0])**2 + (rM[1] - rS[1])**2 + (rM[2] - rS[2])**2
    r = math.sqrt(r2)
    return c/(r**6)

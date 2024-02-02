"""
This script will analyze the WS RMS X,Y sizes as functions of quad gradient
to get transverse Twiss parameters at the entrance of the quad. 

This analysis uses the simplified transport matrices for the quad and a drift
between the quad and the Wire Scanner.

This method based on a parabola cofficients analysis fom our textbook.

At the end, the error analysis is performed 

"""

import sys
import math

from uspas_pylib.matrix_lib import getDriftMatrix2x2

#-------------------------------------------------------------------
#              START of the SCRIPT
#-------------------------------------------------------------------

#---- H- particles mass and energy in MeV
charge = -1.0
mass = 939.9
Ekin = 2.5
momentum = math.sqrt((mass+Ekin)**2 - mass**2)
beta = momentum/(mass+Ekin)
gamma = (mass+Ekin)/mass
momentum /= 1000.
#---- multiplication of magnetic field B and bending radius of curvature Rho
#---- B in Tesla, Rho in meters when momentum is in GeV/c
b_rho = 3.33564*momentum
print ("------------------------------------------------")
print ("Ekin [MeV]        = %7.3f "%Ekin)
print ("Momentum [GeV/c]  = %7.3f "%(momentum))
print ("beta              = ",beta)
print ("gamma             = ",gamma)
print ("B.rho [T*m]      = ",b_rho)
print ("------------------------------------------------")

#---- positions of quad and WS along MEBT in meters
ws_pos = 0.723
quad_pos = 0.418
Lquad = 0.061
#---- this is the distance from the end of the quad to WS
Ldrift = ws_pos - (quad_pos + Lquad/2.0)
print ("------- quad and WS positions ------------------")
print ("quad entrance [m] = %7.3f "%(quad_pos - Lquad/2))
print ("quad center   [m] = %7.3f "%quad_pos)
print ("quad exit     [m] = %7.3f "%(quad_pos + Lquad/2))
print ("WS   center   [m] = %7.3f "%ws_pos)
print ("------------------------------------------------")


#--------------------------------------------------------
#---- Let's read the scan results
#--------------------------------------------------------

#----  sigma2_x_arr = [WS (SigmaX)^2 in m, ...]
#----  sigma2_y_arr = [WS (SigmaY)^2 in m, ...]
#----  K_x_arr = [K0,K1, ... ] - 1/focal_length [1/m]
#----  K_y_arr = [-K0,-K1, ... ] - 1/focal_length [1/m]
sigma2_x_arr = []
sigma2_y_arr = []
K_x_arr = []
K_y_arr = []

with open("rms_data.pkl","rb") as f:
	scan_dic = pickle.load(f)

xrms_arr = scan_dic['xrms']
yrms_arr = scan_dic['yrms']



file_in = open("quad_scan_ws_data.dat","r")
lns = file_in.readlines()
file_in.close()
for ln in lns:
	res_arr = ln.split()
	if(len(res_arr) != 3): continue
	gradient = float(res_arr[0])
	sigmaX = float(res_arr[1])
	sigmaY = float(res_arr[2])
	sigma2_x_arr.append((sigmaX/1000.)**2)
	sigma2_y_arr.append((sigmaY/1000.)**2)
	K = -charge*gradient*Lquad/b_rho
	K_x_arr.append(K)
	K_y_arr.append(-K)
	print (" G [T/m] =  %+7.3f    K [1/m] = %+7.3f   Sigma X,Y = %6.3f , %6.3f"%(gradient,K,sigmaX,sigmaY))

import numpy
import matplotlib.pyplot as plt

#---------------------------------------------
# Let's calculate polynomial fit for X and Y
#---------------------------------------------

#---- polynomial coeff. order reversed in numpy. 
polynom_x, covMtrx_x = numpy.polyfit(K_x_arr,sigma2_x_arr,2,cov=True)

print ("------ X-direction Polynomial coeffs. -----------------------------------")
print (polynom_x)
print ("------ X-direction Variance-Covariant Matrix for Polynomial Coeffs -------")
print (covMtrx_x)

print ("------ X-direction Polynomial coeffs. with errors. ----------------------")
for ind in range(len(polynom_x)):
	print ("x^"+str(ind)," = %+12.5g +-  %10.3g "%(polynom_x[2-ind],math.sqrt(covMtrx_x[2-ind][2-ind])))

poly1d_func = numpy.poly1d(polynom_x)

plt.scatter(K_x_arr,sigma2_x_arr)
plt.plot(K_x_arr, poly1d_func(K_x_arr))
plt.show()

#---- polynomial coeff. order reversed in numpy. 
polynom_y, covMtrx_y = numpy.polyfit(K_y_arr,sigma2_y_arr,2,cov=True)

print ("------ Y-direction Polynomial coeffs. -----------------------------------")
print (polynom_y)
print ("------ Y-direction Variance-Covariant Matrix for Polynomial Coeffs -------")
print (covMtrx_y)

print ("------ Y-direction Polynomial coeffs. with errors. ----------------------")
for ind in range(len(polynom_y)):
	print ("x^"+str(ind)," = %+12.5g +-  %10.3g "%(polynom_y[2-ind],math.sqrt(covMtrx_y[2-ind][2-ind])))

poly1d_func = numpy.poly1d(polynom_y)

plt.scatter(K_y_arr,sigma2_y_arr)
plt.plot(K_y_arr, poly1d_func(K_y_arr))
plt.show()

#---- positions of quad and WS along MEBT in meters
ws_pos = 0.723
quad_pos = 0.418
Lquad = 0.061
#---- this is the distance from the center of the quad to WS
Ldrift = ws_pos - quad_pos

drift_mtrx = getDriftMatrix2x2(Ldrift)
print ("---------------------------------------------------------------------------")
print ("Elements of drift matrix:")
print ("S11 = %+12.5g "%drift_mtrx.get(0,0))
print ("S12 = %+12.5g "%drift_mtrx.get(0,1))
print ("S22 = %+12.5g "%drift_mtrx.get(1,1))
print ("---------------------------------------------------------------------------")

#----------------------------------
# To be contunued
#----------------------------------
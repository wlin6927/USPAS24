"""
This script will analyze the WS RMS X,Y sizes as functions of quad gradient
to get transverse Twiss parameters at the entrance of the quad. 

This analysis uses the simplified transport matrices for the quad and a drift
between the quad and the Wire Scanner.

There are two possible representation of a quad transport matrix
1. drift-thin_quad-drift
2. non-zero length quad

At the end, the error analysis is performed using LSQ methods.

"""

import sys
import math
import random
import time
import pickle

from orbit.core.orbit_utils import Matrix
from orbit.core.orbit_utils import PhaseVector

from uspas_pylib.matrix_lib import printMatrix
from uspas_pylib.matrix_lib import printVector
from uspas_pylib.matrix_lib import getDriftMatrix2x2
from uspas_pylib.matrix_lib import getThinQuadMatrix2x2
from uspas_pylib.matrix_lib import getQuadMatrix2x2
from uspas_pylib.matrix_lib import getCovarianceMatrix
from uspas_pylib.matrix_lib import getCovarianceMatrixWithWeights
from uspas_pylib.matrix_lib import getBeamCorrelations
from uspas_pylib.matrix_lib import getBeamCorrelationsWithWeights
from uspas_pylib.matrix_lib import getErrorWeghtMatrix

from uspas_pylib.twiss_parameters_lib import emittanceFunc
from uspas_pylib.twiss_parameters_lib import betaFunc
from uspas_pylib.twiss_parameters_lib import alphaFunc

def getTransportMatrix(momentum,G,Lquad,Ldrift,charge, quadMatrixGenerator = getQuadMatrix2x2):
	"""
	Returns PyORBIT matrix which is a transport matrix for (x_new,xp_new)^T = M * (x,xp)^T
	for "quad-drift" lattice.
	momentum - in GeV/c
	G - quad gradient in T/m
	Lquad - quad legth in m
	Ldrift - drift length in m
	charge - +1 (proton) or -1 (H- ion)
	quadMatrixGenerator could be getQuadMatrix2x2 or getThinQuadMatrix2x2
	"""
	drift_mtr = getDriftMatrix2x2(Ldrift)
	q_mtr = quadMatrixGenerator(momentum,G,Lquad,charge)
	#---- Lquad != 0. - thick quad transport matrix
	matr = drift_mtr.mult(q_mtr)
	return matr

def getLSQ_Matrix(gradient_arr,momentum,Lquad,Ldrift,charge, quadMatrixGenerator = getQuadMatrix2x2):
	"""
	Creates the LSQ matrix from the 1st lines of transport matrix.
	Matrix defines the system of equations:
	RMS^2[ind] = mtrx[ind][0]*<x^2> + mtrx[ind][1]*<x*x'> +  mtrx[ind][2]*<x'^2> 
	where ind = 0...(len(results_arr)-1)
	<x^2> , <x*x'>, <x'^2> - unknown correlations we want to find. 
	"""
	n_row = len(gradient_arr)
	mtrx = Matrix(n_row,3)
	for ind in range(n_row):
		G = gradient_arr[ind]
		transport_matrix = getTransportMatrix(momentum,G,Lquad,Ldrift,charge, quadMatrixGenerator)
		m1 = transport_matrix.get(0,0)
		m2 = transport_matrix.get(0,1)
		mtrx.set(ind,0,m1**2)
		mtrx.set(ind,1,2*m1*m2)
		mtrx.set(ind,2,m2**2)
	return mtrx

def getRMS2_Vector(results_arr):
	"""
	Returns PyORBIT PhaseVector with RMS^2 sizes from Wire Scanner
	for different quadrupole gradients.
	results_arr = [[G,rms^2], ...]
	"""
	n_row = len(results_arr)
	rms2_vector = PhaseVector(n_row)
	for ind, [G, rms] in enumerate(results_arr):
		rms2_vector.set(ind,rms**2)
	return rms2_vector
 
def getError(func,param_val_vector,covariance_mtrx, step_arr = [0.01,0.01,0.01]):
	"""
	Calculates error for function of 3 variables using known values of these variables
	and known variance-covariance matrix for errors of variables.
	"""
	vct0 = param_val_vector
	derivative_arr = []
	for ind in range(3):
		delta = abs(vct0.get(ind))*step_arr[ind]
		vct1 = PhaseVector(param_val_vector)
		vct1.set(ind,vct1.get(ind) + delta)
		derivative = (func(vct1) - func(vct0))/delta
		derivative_arr.append(derivative)
	error2 = 0.
	for ind0 in range(3):
		for ind1 in range(3):
			error2 += derivative_arr[ind0]*derivative_arr[ind1]*covariance_mtrx.get(ind0,ind1)
	return math.sqrt(error2)
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

#----  results_x_arr = [[G in T/m, WS SigmaX in mm], ...]
#----  results_y_arr = [[G in T/m, WS SigmaY in mm], ...]
#----  G_x_arr = [G0,G1, ... ] - G in T/m
#----  G_y_arr = [-G0,-G1, ... ] - G in T/m
results_x_arr = []
results_y_arr = []
G_x_arr = []
G_y_arr = []

with open("rms_data1.pkl","rb") as f:
	scan_dic = pickle.load(f)

G_x_arr = scan_dic['QH03']
G_y_arr = -1 * scan_dic['QH03']

xrms_arr = scan_dic['xrms']
yrms_arr = scan_dic['yrms']

for i in range(len(G_x_arr)):
	G = G_x_arr[i]
	sigmaX = xrms_arr[i]
	sigmaY = yrms_arr[i]
	results_x_arr.append([G, sigmaX])
	results_y_arr.append([G, sigmaY])
	print(" G [T/m] = %+7.3f   Sigma X,Y = %6.3f , %6.3f" % (G, sigmaX, sigmaY))


rms2_x_vector = getRMS2_Vector(results_x_arr)
rms2_y_vector = getRMS2_Vector(results_y_arr)

#--------------------------------------------------------------------
#   Calculations of [<x*x>,<x*x'>,<x'*x'>] or [<y*y>,<y*y'>,<y'*y'>]
#-------------------------------------------------------------------

quadMatrixGenerator = getQuadMatrix2x2
#quadMatrixGenerator = getThinQuadMatrix2x2

lsq_x_matrix = getLSQ_Matrix(G_x_arr,momentum,Lquad,Ldrift,charge,quadMatrixGenerator)
lsq_y_matrix = getLSQ_Matrix(G_y_arr,momentum,Lquad,Ldrift,charge,quadMatrixGenerator)

#---- V(results) = (M^T * M)^-1 * M^T * V(rms2) where M is LSQ_Matrix
#---- Results are [<x*x>,<x*x'>,<x'*x'>] and [<y*y>,<y*y'>,<y'*y'>]
#---- in the form of PyORBIT PhaseVerctors

print("========Horizontal vector from LSQ=======")
lsq_x_matrix_T = lsq_x_matrix.copy()
lsq_x_matrix_T.transpose()

big_matrix = lsq_x_matrix_T.copy()
result_matrix = big_matrix.mult(lsq_x_matrix)

invert_matrix = result_matrix.invert()
inv_trans_matrix = invert_matrix.mult(lsq_x_matrix_T)
corr_x_vector = inv_trans_matrix.mult(rms2_x_vector)
printVector(corr_x_vector)

alpha_x = -corr_x_vector.get(1)/math.sqrt(corr_x_vector.get(0)*corr_x_vector.get(2) - (corr_x_vector.get(1))**2)
print("alpha_x:",alphaFunc(corr_x_vector))
beta_x = corr_x_vector.get(0)/math.sqrt(corr_x_vector.get(0)*corr_x_vector.get(2) - (corr_x_vector.get(1))**2)
print("beta_x:",betaFunc(corr_x_vector))
emit_x = math.sqrt(corr_x_vector.get(0)*corr_x_vector.get(2) - (corr_x_vector.get(1))**2)
print("emit_x:",emittanceFunc(corr_x_vector))

print("========Vertical vector from LSQ=======")
lsq_y_matrix_T = lsq_y_matrix.copy()
lsq_y_matrix_T.transpose()
big_matrix_y = lsq_y_matrix_T.copy()
result_matrix_y = big_matrix_y.mult(lsq_y_matrix)

invert_matrix_y = result_matrix_y.invert()
inv_trans_matrix_y = invert_matrix_y.mult(lsq_y_matrix_T)
corr_y_vector = inv_trans_matrix_y.mult(rms2_y_vector)
printVector(corr_y_vector)

alpha_y = -corr_y_vector.get(1)/math.sqrt(corr_y_vector.get(0)*corr_y_vector.get(2) - (corr_y_vector.get(1))**2)
print("alpha_y:",alphaFunc(corr_y_vector))
beta_y = corr_y_vector.get(0)/math.sqrt(corr_y_vector.get(0)*corr_y_vector.get(2) - (corr_y_vector.get(1))**2)
print("beta_y:",betaFunc(corr_y_vector))
emit_y = math.sqrt(corr_y_vector.get(0)*corr_y_vector.get(2) - (corr_y_vector.get(1))**2)
print("emit_y:",emittanceFunc(corr_y_vector))

#-----------------------------------------
# To be continued
#-----------------------------------------

sys.exit(0)


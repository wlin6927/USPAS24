"""
This script is the start of the home work where you
have to close 3-kickers bump using DCV01, DCV04, 
and DCV05 correctors in MEBT.

It shows how to create PV names and communicate with 
the Virtual Accelerator.

>virtual_accelerator --debug  --sequences MEBT

"""

import sys
import math
import time
import numpy as np

from epics import pv as pv_channel

# pyorbit optimization packages
from orbit.utils.fitting import Solver
from orbit.utils.fitting import Scorer
from orbit.utils.fitting import SolveStopperFactory
from orbit.utils.fitting import ScoreboardActionListener
from orbit.utils.fitting import VariableProxy
from orbit.utils.fitting import TrialPoint
from orbit.utils.fitting import SimplexSearchAlgorithm

#-------------------------------------------------------------------
#              START of the SCRIPT
#-------------------------------------------------------------------

#---- DCV01, DCV04, DCV05, DCV10, DCV11, DCV14  
dcv_ind_arr = [1,4,5]
dcv_field_pv_arr = []
for ind in dcv_ind_arr:
	dcv_field_pv_arr.append(pv_channel.PV("MEBT_Mag:PS_DCV"+"%02d"%ind+":B_Set"))

#---- put 0. [T] field in all DCV
for pv in dcv_field_pv_arr:
	pv.put(0.)

time.sleep(1.0)

#---- BPM10, BPM11, BPM14
bpm_ind_arr = [1,4,5,10,11,14]
bpm_ver_pos_pv_arr = []
for ind in bpm_ind_arr:
	bpm_ver_pos_pv_arr.append(pv_channel.PV("MEBT_Diag:BPM"+"%02d"%ind+":yAvg"))

#-------------------------------------------------
#---- Let's print BPM signals before bump closing
#-------------------------------------------------
print ("======= Vertical orbit [mm] before Closing =====")
bpm_orbit = []
for bpm_pv in bpm_ver_pos_pv_arr:
	bpm_orbit.append(bpm_pv.get())
print (bpm_orbit)

class OrbitFittingScorer(Scorer):
	"""
	The Scorer implementation for closing bump orbit
	"""

	def __init__(self, dcv04_field, dcv05_field):
		self.dcv04_field = dcv04_field
		self.dcv05_field = dcv05_field


	def getScore(self, trialPoint):
		"""
		Implementation of getScore method for Scorer class.
		Calculates sum of squares of differences between goal and fitted orbit values
		"""
		[dcv04_field, dcv05_field] = trialPoint.getVariableProxyValuesArr()
		dcv_field_pv_arr[1].put(dcv04_field)
		dcv_field_pv_arr[2].put(dcv05_field)
		time.sleep(2.0)
		diff2 = 0
		for bpm_pv in bpm_ver_pos_pv_arr[3:]:
			bpm_orbit = bpm_pv.get()
			diff2 += (bpm_orbit - 0.) ** 2
		return diff2

def orbit_fit(dcv01_field, dcv04_field, dcv05_field, dcv_step, maxIter):

	# ---- set DCV01 t0 field
	print("Initial DCV01 field set to",dcv01_field,"Tesla.")
	dcv_field_pv_arr[0].put(dcv01_field)

	# ---- give the accelerator time to put beam through MEBT
	time.sleep(2.0)

	print("======= Vertical orbit [mm] before Closing =====")
	bpm_orbit = []
	for bpm_pv in bpm_ver_pos_pv_arr:
		bpm_orbit.append(bpm_pv.get())
	print(bpm_orbit)

	# Trial point
	trialPoint = TrialPoint()
	variableProxy = VariableProxy("DCV04", dcv04_field, dcv_step)
	trialPoint.addVariableProxy(variableProxy)
	variableProxy = VariableProxy("DCV05", dcv05_field, dcv_step)
	trialPoint.addVariableProxy(variableProxy)

	# Instance of Scorer class
	scorer = OrbitFittingScorer(dcv04_field,dcv05_field)

	# Search algorithm from PyORBIT native package
	searchAlgorithm = SimplexSearchAlgorithm()

	solverStopper = SolveStopperFactory.maxIterationStopper(maxIter)

	solver = Solver()
	solver.setAlgorithm(searchAlgorithm)
	solver.setStopper(solverStopper)

	print("======= START orbit closing =======")
	# Fitting process itself
	solver.solve(scorer, trialPoint)

	trialPoint = solver.getScoreboard().getBestTrialPoint()
	print (trialPoint.textDesciption())
	print ("=======Orbit closing END==========")

	print("======= Vertical orbit [mm] after Closing =====")
	bpm_orbit = []
	for bpm_pv in bpm_ver_pos_pv_arr:
		bpm_orbit.append(bpm_pv.get())
	print(bpm_orbit)

	dcv04_field = pv_channel.PV("MEBT_Mag:PS_DCV04:B").get()
	dcv05_field = pv_channel.PV("MEBT_Mag:PS_DCV05:B").get()

	# effective length of correctors
	dcv01_len = 0.061
	dcv04_len = 0.061
	dcv05_len = 0.066
	print("======= Relations between correctors’ kicks =======")
	print("kick1/kick2:", dcv04_field * dcv04_len / (dcv01_field * dcv01_len))
	print("kick1/kick3:", dcv05_field * dcv05_len / (dcv01_field * dcv01_len))

	print("======= Restoring all correctors to zero =======")
	# ---- put 0. [T] field in all DCV
	for pv in dcv_field_pv_arr:
		pv.put(0.)
	time.sleep(2.0)
	return

dcv_step = 0.001 # Tesla
dcv04_field = pv_channel.PV("MEBT_Mag:PS_DCV04:B").get()
dcv05_field = pv_channel.PV("MEBT_Mag:PS_DCV05:B").get()


maxIter = 70

dcv01_fields = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08]

for dcv01_field in dcv01_fields:
	orbit_fit(dcv01_field, dcv04_field, dcv05_field, dcv_step, maxIter)


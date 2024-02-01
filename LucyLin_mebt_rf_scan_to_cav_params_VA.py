"""
This script is the start of the home work where you
will scan a phase of the last MEBT re-buncher RF caity
to find the non-accelerating phase and the cavity amplitude.

>virtual_accelerator --debug  --sequences MEBT

"""

import os
import sys
import math
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from epics import pv as pv_channel

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory
from orbit.py_linac.lattice import MarkerLinacNode
from orbit.utils import phaseNearTargetPhaseDeg

from orbit.core.orbit_utils import Function
from orbit.utils.fitting import PolynomialFit

from uspas_pylib.harmonic_data_fitting_lib import fitCosineFunc

#-------------------------------------------------------------------
#              START of the SCRIPT
#-------------------------------------------------------------------

#-----Parameters at the entrance of MEBT ---------------
# transverse emittances are unnormalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 0.0025 # in [GeV]
mass = 0.938272    # in [GeV]
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
print ("relat. gamma=",gamma)
print ("relat.  beta=",beta)
bpm_frequency = 805.0e+6 # MHz
v_light = 2.99792458e+8  # in [m/sec]



names = ["MEBT",]
#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()

#---- the XML file name with the structure
xml_file_name = os.environ["HOME"] + "/uspas24-CR/lattice/sns_linac.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

#---- dictionary with the start and end positions of 1st level nodes
node_position_dict = accLattice.getNodePositionsDict()

rf_cavs = accLattice.getRF_Cavities()

#---- we will collect all BPMs accelerator nodes even they are child nodes
bpms = []
for node in accLattice.getNodes():
	if(isinstance(node,MarkerLinacNode) and node.getName().find("BPM") >= 0):
		bpms.append(node)
		continue
	for childNode in node.getBodyChildren():
		if(isinstance(childNode,MarkerLinacNode) and childNode.getName().find("BPM") >= 0):
			bpms.append(childNode)
			continue

#---- now we get the positions of the last cavity and the last BPM
cav_index = len(rf_cavs) - 1
bpm_index = len(bpms) - 1
cav = rf_cavs[cav_index]
bpm = bpms[bpm_index]
bpm_pos = bpm.getPosition()
cav_pos = cav.getPosition()
L_dist = bpm_pos - cav_pos 
print ("========================================")
print ("Cavity = %15s "%cav.getName()," pos[m] = %6.3f "%cav_pos)
print ("BPM    = %15s "%bpm.getName()," pos[m] = %6.3f "%bpm_pos)
print ("========================================")

#------------------------------------------------------------------
#---- Now we perform a phase scan measuring BPM phase
#------------------------------------------------------------------
cav_phase_pv = pv_channel.PV("MEBT_LLRF:FCM" + str(cav_index+1) + ":CtlPhaseSet")
cav_amp_pv   = pv_channel.PV("MEBT_LLRF:FCM" + str(cav_index+1) + ":CtlAmpSet")
bpm_phase_pv = pv_channel.PV("MEBT_Diag:BPM14:phaseAvg")
bpm_amp_pv   = pv_channel.PV("MEBT_Diag:BPM14:amplitudeAvg")

bpm_phase0 = bpm_phase_pv.get()
cav_phase0 = cav_phase_pv.get()

print ("========================================")
print ("Cavity = %15s "%cav.getName()," phase[deg] = %6.3f "%cav_phase0)
print ("BPM    = %15s "%bpm.getName()," phase[deg] = %6.3f "%bpm_phase0)
print ("========================================")

cav_phases = np.linspace(cav_phase0 - 180, cav_phase0 + 180, 12)

bpm_phases = []
for phase in cav_phases:
    cav_phase_pv.put(phase)
    time.sleep(1.0)
    bpm_phase = bpm_phase_pv.get()
    bpm_phases.append(bpm_phase)

# return to original phase
cav_phase_pv.put(cav_phase0)

# plt.plot(cav_phases,bpm_phases,'--o')
# plt.xlabel('Cavity phase [deg]')
# plt.ylabel('BPM phase [deg]')
# plt.show()

print("============uspas_pylib fitCosineFunc fitting============")
((amp,phase_offset,avg_val),scorer) = fitCosineFunc(cav_phases, bpm_phases, show_progress = False)
print ("(phase_offset,amp,avg_val) = ",(phase_offset,amp,avg_val))
print ("============Stop.============")

dt = amp / 360 / bpm_frequency
t0 = L_dist / (beta * v_light)

v = L_dist / (t0 - dt)
beta1 = v / v_light
gamma1 = 1 / (math.sqrt(1 - beta1*beta1))


e_kin_new = (gamma1 - 1) * mass

Veff = (e_kin_new - e_kin_ini) * 1e9
print("V_eff of the cavity from 360 phase scan is:", Veff/1000, "kV.")


cav_phases1 = np.linspace(91.7, 93.5, 7)
cav_amps = np.linspace(0.2,1,5)

bpm_phases = []

for amp in cav_amps:
    cav_amp_pv.put(amp)
    time.sleep(1)
    bpm_phases1 = []
    for phase in cav_phases1:
        cav_phase_pv.put(phase)
        time.sleep(1)
        bpm_phase = bpm_phase_pv.get()
        bpm_phases1.append(bpm_phase)
    bpm_phases.append(bpm_phases1)

# plt.plot(cav_phases1,bpm_phases[0],'--o')
# plt.plot(cav_phases1,bpm_phases[1],'--o')
# plt.plot(cav_phases1,bpm_phases[2],'--o')
# plt.plot(cav_phases1,bpm_phases[3],'--o')
# plt.plot(cav_phases1,bpm_phases[4],'--o')
# plt.axvline(x=92.19,color='k')
# plt.xlabel("Cavity phase [deg]")
# plt.ylabel("BPM phase [deg]")
# plt.savefig("zero_crossing.png",bbox_inches ="tight");

slope = np.mean(np.diff(bpm_phases[4]) / np.diff(cav_phases1))
print("Slope of dphi_BPM / dphi_RF:", slope)

Veff = - slope * v_light * beta * 2 * e_kin_ini*1e9 / (2*np.pi*bpm_frequency * L_dist)
print("V_eff of the cavity from zero-crossing scan is:", Veff/1000, "kV.")



sys.exit(0)
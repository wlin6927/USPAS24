"""
This script will create the MEBT model and will synchronize it with VA.
Then, it will calculate transport matrices between DCV01, DCV04, and DCV05,
and it will calulate the fields in  DCV04 and DCV05 to close the bump
created by DCV01.

Finally, it will check the results with VA.

>virtual_accelerator --debug  --sequences MEBT

"""

import os
import sys
import math
import time
import numpy as np
import matplotlib.pyplot as plt

from epics import pv as pv_channel

from orbit.py_linac.linac_parsers import SNS_LinacLatticeFactory
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.bunch_generators import TwissContainer

from orbit.lattice import AccNode
from orbit.py_linac.lattice import MarkerLinacNode
from orbit.py_linac.lattice import DCorrectorH, DCorrectorV

from orbit.py_linac.lattice import LinacTrMatricesContrioller

from orbit.core.orbit_utils import Matrix

from uspas_pylib.sns_linac_bunch_generator import SNS_Linac_BunchGenerator
from uspas_pylib.bpm_model_node_lib import ModelBPM
from uspas_pylib.bpm_model_node_lib import addModelBPM

from uspas_pylib.matrix_lib import printMatrix
from uspas_pylib.matrix_lib import printVector

#---- instance of a class for Bunch analysis
twiss_analysis = BunchTwissAnalysis()

#-------------------------------------------------------------------
#              START of the SCRIPT
#-------------------------------------------------------------------

#--------------------------------------------
#---- Let's make the lattice
#--------------------------------------------

names = ["MEBT",]
#---- create the factory instance
sns_linac_factory = SNS_LinacLatticeFactory()

#---- the XML file name with the structure
xml_file_name = os.environ["HOME"] + "/uspas24-CR/lattice/sns_linac.xml"

#---- make lattice from XML file 
accLattice = sns_linac_factory.getLinacAccLattice(names,xml_file_name)

#---- dictionary with the start and end positions of 1st level nodes
node_position_dict = accLattice.getNodePositionsDict()

#---------------------------------------------
#---- Let's collect BPM, DCV, and Quad nodes
#---------------------------------------------

bpm_model_nodes = []
for node in accLattice.getNodes():
    position = node.getPosition()
    if(isinstance(node,MarkerLinacNode) and node.getName().find("BPM") >= 0):
        bpm_model_node = addModelBPM(node,position)
        bpm_model_nodes.append(bpm_model_node)
        #print ("debug bpm-model=",bpm_model_node.getName()," position=",bpm_model_node.getPosition())
        continue
    for childNode in node.getBodyChildren():
        if(isinstance(childNode,MarkerLinacNode) and childNode.getName().find("BPM") >= 0):
            bpm_model_node = addModelBPM(childNode,position)
            bpm_model_nodes.append(bpm_model_node)
            #print ("debug bpm-model=",bpm_model_node.getName()," position=",bpm_model_node.getPosition())
            continue

dcv_nodes = []
for node in accLattice.getNodes():
    position = node.getPosition()
    if(isinstance(node,DCorrectorV)):
        dcv_nodes.append(node)
        #print ("debug dvc=",dcv_node.getName()," position=",dcv_node.getPosition())
        continue
    for childNode in node.getBodyChildren():
        if(isinstance(childNode,DCorrectorV)):
            dcv_nodes.append(childNode)
            #print ("debug dcv=",childNode.getName()," position=",childNode.getPosition())
            continue

quad_nodes = accLattice.getQuads()

#----------------------------------------------------------------------
#---- Now let's synchronize PyORBIT model with Virtual Accelerator 
#----------------------------------------------------------------------
print ("====================================")
for quad_node in quad_nodes:
    pv = pv_channel.PV(quad_node.getName().replace(":",":PS_")+":B_Set")
    field_grad = pv.get()
    quad_node.setField(field_grad)
    print ("debug quad=",quad_node.getName()," G[T/m] =",field_grad)
print ("====================================")
dcv_field_pv_arr = []
for dcv_node in dcv_nodes:
    pv = pv_channel.PV(dcv_node.getName().replace(":",":PS_")+":B_Set")
    dcv_field_pv_arr.append(pv)
    pv.put(0.)
    print ("debug dcv=",dcv_node.getName()," pv=",pv.pvname)
print ("====================================")
bpm_ver_pos_pv_arr = []
for bpm_model_node in bpm_model_nodes:
    pv = pv_channel.PV(bpm_model_node.getName().replace("-model","")+":yAvg")
    bpm_ver_pos_pv_arr.append(pv)
    yAvg = pv.get()
    print ("debug bpm pv =",pv.pvname," yAvg[mm]= %+6.4f"%yAvg)
print ("====================================")

#----------------------------------------------------------------------
#---- Now we generate the bunch for the MEBT model at the entrance 
#----------------------------------------------------------------------

#-----TWISS Parameters at the entrance of MEBT ---------------
# transverse emittances are unnormalized and in pi*mm*mrad
# longitudinal emittance is in pi*eV*sec
e_kin_ini = 0.0025 # in [GeV]

#------ PyORBIT MEBT entrance Twiss parameters 
(alphaX,betaX,emittX) = (-1.9620,   0.1831, 2.8764e-6)
(alphaY,betaY,emittY) = ( 1.7681,   0.1620, 2.8764e-6)
(alphaZ,betaZ,emittZ) = (-0.0196, 116.4148, 0.0165e-6)
 
twissX = TwissContainer(alphaX,betaX,emittX)
twissY = TwissContainer(alphaY,betaY,emittY)
twissZ = TwissContainer(alphaZ,betaZ,emittZ)

bunch_gen = SNS_Linac_BunchGenerator(twissX,twissY,twissZ)
#set the initial kinetic energy in GeV
bunch_gen.setKinEnergy(e_kin_ini)

bunch_in = bunch_gen.getBunch(nParticles = 10000)
#bunch_in = bunch_gen.getBunch(nParticles = 10000)
#bunch_in = bunch_gen.getBunch(nParticles = 10000)
print ("Bunch is ready. N particles=",bunch_in.getSize())
print ("====================================")

#--------------------------------------------------------------------
#---- Setup the nodes for transport matrix calculations
#--------------------------------------------------------------------
trMatricesGenerator = LinacTrMatricesContrioller()
#---- from DCV01, DCV04, DCV05, DCV10, DCV11, DCV14 we use 1st 3 correctors
lattice_trMatrix_nodes = [dcv_nodes[0],dcv_nodes[1],dcv_nodes[2]]
#---- these matrix nodes will keep transport matrices after bunch tracking
trMatrixNodes = trMatricesGenerator.addTrMatrxGenNodesAtEntrance(accLattice,lattice_trMatrix_nodes)

#---- The use of Twiss weights makes transport matrices more accurate.
for trMtrxNode in trMatrixNodes:
    #---- setting up to use Twiss weights for transport matrix calculations
    #---- for X,Y,Z directions.
    trMtrxNode.setTwissWeightUse(True,True,True)

#----------------------------------------------------------------------
#---- Now we track bunch through MEBT to get transport matrices
#----------------------------------------------------------------------

#---- we will keep bunch_in and always start with its copy
bunch = Bunch()
bunch_in.copyBunchTo(bunch)

#---- set up design for RF cavities
accLattice.trackDesignBunch(bunch)

#---- really track the bunch through MEBT
accLattice.trackBunch(bunch)

#---------------------------------------------------------------------
#---- Now we print transport matrices between DCV01, DCV04, DCV05
#---------------------------------------------------------------------
#---- There are 3 matrix 
#---- 1. from DCV01 to DCV01 - unit matrix
#---- 2. from DCV01 to DCV04
#---- 3. from DCV01 to DCV05
for trMtrxNode in trMatrixNodes:
    print ("Transport Matrix Node =",trMtrxNode.getName()," pos[m] = %5.3f"%trMtrxNode.getPosition())
    trMtrx = trMtrxNode.getTransportMatrix()
    printMatrix(trMtrx)
    print ("====================================")

#---- here our transport matrices between correctors 
trMtrx_1_to_4 = trMatrixNodes[1].getTransportMatrix()
trMtrx_1_to_5 = trMatrixNodes[2].getTransportMatrix()
trMtrx_4_to_5 = trMtrx_1_to_5.mult(trMtrx_1_to_4.invert())

for dcv_node in dcv_nodes:
    print ("Corrector =",dcv_node.getName()," L[m]=",dcv_node.getParam("effLength"))
print ("====================================")

#---- length of correctors
length_1 = dcv_nodes[0].getParam("effLength")
length_4 = dcv_nodes[1].getParam("effLength")
length_5 = dcv_nodes[2].getParam("effLength")

#---------------------------------------------------------------
#---- Now we change the DCV01 field to 0.01 T and see the BPMs
#---------------------------------------------------------------
dc01_field = 0.01
dcv_field_pv_arr[0].put(dc01_field)

time.sleep(2.1)

for pv in bpm_ver_pos_pv_arr:
    yAvg = pv.get()
    print ("debug bpm pv =",pv.pvname," yAvg[mm]= %+6.4f"%(yAvg))
print ("====================================")

# calculate beam profile at DCV04
m12_1to4 = trMtrx_1_to_4.get(2,3)
m22_1to4 = trMtrx_1_to_4.get(3,3)

m11_4to5 = trMtrx_4_to_5.get(2,2)
m12_4to5 = trMtrx_4_to_5.get(2,3)
m21_4to5 = trMtrx_4_to_5.get(3,2)
m22_4to5 = trMtrx_4_to_5.get(3,3)

dc04_field = -(m11_4to5 * m12_1to4 + m12_4to5 * m22_1to4)\
             * dc01_field * length_1 / (m12_4to5 * length_4)
print("DC04 field:", dc04_field)

dc05_field = -1/length_5 * ((m21_4to5*m12_1to4 + m22_1to4*m22_4to5)\
                            *dc01_field*length_1+m22_4to5*dc04_field*length_4)
print("DC05 field:", dc05_field)

dcv_field_pv_arr[1].put(dc04_field)
dcv_field_pv_arr[2].put(dc05_field)

time.sleep(2.1)

for pv in bpm_ver_pos_pv_arr:
    yAvg = pv.get()
    print ("bpm pv =",pv.pvname," yAvg[mm]= %+6.4f"%(yAvg))
print ("====================================")

print("Now track beam along the line")

#-------------------------------------------------------------
# How to see bunch parameters along the tracking
#-------------------------------------------------------------
from orbit.lattice import AccActionsContainer

# track through the lattice
paramsDict = {"old_pos": -1.0, "count": 0, "pos_step": 0.1}
actionContainer = AccActionsContainer("Bunch Tracking")

pos_start = 0.0

node_positions = []
node_y = []

def action_entrance(paramsDict):
    node = paramsDict["node"]
    bunch = paramsDict["bunch"]
    pos = paramsDict["path_length"]
    if paramsDict["old_pos"] == pos:
        return
    if paramsDict["old_pos"] + paramsDict["pos_step"] > pos:
        return
    paramsDict["old_pos"] = pos
    paramsDict["count"] += 1
    twiss_analysis.analyzeBunch(bunch)
    y = twiss_analysis.getAverage(2) * 1000 # m to mm

    node_positions.append(pos)
    node_y.append(y)

def action_exit(paramsDict):
    action_entrance(paramsDict)


actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

bunch = Bunch()
bunch_in.copyBunchTo(bunch)

# set DCV01 to 0.01 T
for node in dcv_nodes:
    if "DCV01" in node.getName():
        node.setField(dc01_field)
    else:
        node.setField(0.0)

accLattice.trackBunch(bunch, paramsDict = paramsDict, actionContainer = actionContainer)

pos_bef = np.array(node_positions)
y_bef = np.array(node_y)

paramsDict = {"old_pos": -1.0, "count": 0, "pos_step": 0.1}
actionContainer = AccActionsContainer("Bunch Tracking")

pos_start = 0.0

node_positions = []
node_y = []

pos_start = 0.0

actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)
actionContainer.addAction(action_exit, AccActionsContainer.EXIT)

bunch2 = Bunch()
bunch_in.copyBunchTo(bunch2)


# set DCV01 to 0.01 T
for node in dcv_nodes:
    name = node.getName()
    if "DCV01" in name:
        node.setField(dc01_field)
    elif "DCV04" in name:
        node.setField(dc04_field)
    elif "DCV05" in name:
        node.setField(dc05_field)
    else:
        node.setField(0.0)


accLattice.trackBunch(bunch2, paramsDict = paramsDict, actionContainer = actionContainer)

pos_opt = np.array(node_positions)
y_opt = np.array(node_y)

import pickle
dic = {}
dic['pos_bef'] = pos_bef
dic['y_bef'] = y_bef
dic['pos_opt'] = pos_opt
dic['y_opt'] = y_opt

with open("data.pkl","wb") as f:
    pickle.dump(dic,f)

plt.plot(pos_bef,y_bef,'--o',label='before')
plt.plot(pos_opt,y_opt,'--o',label='after')
plt.xlabel("position [m]")
plt.ylabel("Vertical orbit [mm]")
plt.legend()
plt.show()
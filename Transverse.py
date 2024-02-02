import os
import sys
import math
import random
import time

import pickle

import numpy as np
import matplotlib.pyplot as plt

from epics import pv as pv_channel

# positions of quad and WS centers along MEBT in meters
ws_pos = 0.723
quad_pos = 0.418
Lquad = 0.061

qv04_pv = pv_channel.PV("MEBT_Mag:PS_QV04:B_Set") # should be set to 0.
rf1_pv = pv_channel.PV("MEBT_LLRF:FCM1:CtlAmpSet") # should be set to 0.
qh03_pv = pv_channel.PV("MEBT_Mag:PS_QH03:B_Set") # quad gradient in the range 4-40 T/m
ws_pos_pv = pv_channel.PV("MEBT_Diag:WS04a:Position")
ws_pos_set_pv = pv_channel.PV("MEBT_Diag:WS04a:Position_Set")
ws_speed_set_pv = pv_channel.PV("MEBT_Diag:WS01:Speed_Set")
ws_hcont_pv = pv_channel.PV("MEBT_Diag:WS04a:Hor_Cont")
ws_vcont_pv = pv_channel.PV("MEBT_Diag:WS04a:Ver_Cont")

qh03_range = np.linspace(-40,-4,10)

# set qv04 and rf1 to 0 to make drift section
qv04_pv.put(0.0)
rf1_pv.put(0.0)
time.sleep(1.0)

# qh03_pv.put(qh03)

print(qv04_pv.get(),rf1_pv.get())
ws_speed_set_pv.put(3.0)
time.sleep(1.0)

vcont_master_arr = []
hcont_master_arr = []

# print("Start vertical wire scan:")
# ws_pos_vrange = np.linspace(-25,-5,40)
# vcont_arr = []
# for ws_pos in ws_pos_vrange:
#     ws_pos_set_pv.put(ws_pos)
#     while abs(ws_pos_pv.get() - ws_pos) > 0.1:
#         time.sleep(1.0)
#     time.sleep(1.5)
#     vcont = ws_vcont_pv.get()
#     vcont_arr.append(vcont)
#
# print("Start horizontal wire scan:")
# ws_pos_hrange = np.linspace(10,20,40)
# hcont_arr = []
# for ws_pos in ws_pos_hrange:
#     ws_pos_set_pv.put(ws_pos)
#     while abs(ws_pos_pv.get() - ws_pos) > 0.1:
#         time.sleep(1.0)
#     time.sleep(1.5)
#     hcont = ws_hcont_pv.get()
#     hcont_arr.append(hcont)

# plt.plot(ws_pos_vrange,vcont_arr)
# plt.plot(ws_pos_hrange,hcont_arr)
# plt.show()

count = 0

for qh03 in qh03_range:

    print("Put wire scanner back...")
    ws_speed_set_pv.put(100)
    ws_pos_set_pv.put(-25)
    time.sleep(2.0)
    ws_speed_set_pv.put(3.0)
    time.sleep(1.0)

    print("QH03 set to field:",qh03)
    qh03_pv.put(qh03)
    time.sleep(1.5)


    print("Start vertical wire scan for QH03 field:", qh03)
    ws_pos_vrange = np.linspace(-25,-5,40)
    vcont_arr = []
    for ws_pos in ws_pos_vrange:
        ws_pos_set_pv.put(ws_pos)
        while abs(ws_pos_pv.get() - ws_pos) > 0.1:
            time.sleep(1.0)
        time.sleep(1.5)
        vcont = ws_vcont_pv.get()
        vcont_arr.append(vcont)

    vcont_master_arr.append(vcont_arr)

    print("Start horizontal wire scan for QH03 field:", qh03)
    ws_pos_hrange = np.linspace(10,20,40)
    hcont_arr = []
    for ws_pos in ws_pos_hrange:
        ws_pos_set_pv.put(ws_pos)
        while abs(ws_pos_pv.get() - ws_pos) > 0.1:
            time.sleep(1.0)
        time.sleep(1.5)
        hcont = ws_hcont_pv.get()
        hcont_arr.append(hcont)

    hcont_master_arr.append(hcont_arr)

    minidic = {}
    minidic['QH03'] = qh03
    minidic['Hcont'] = hcont_arr
    minidic['Vcont'] = vcont_arr

    with open("scan_data1_"+str(count)+".pkl", "wb") as f:
        pickle.dump(minidic, f)

    count += 1

# dic = {}
# dic['QH03'] = qh03_range
# dic['Hcont'] = hcont_master_arr
# dic['Vcont'] = vcont_master_arr
#
# with open("scan_data.pkl","wb") as f:
#     pickle.dump(dic,f)








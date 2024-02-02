"""
This script will perform the Wire Scans WS04a for different gradients
of the MEBT:QH03.

It shows how to create PV names and communicate with 
the Virtual Accelerator.

Start of Virtual Accelerator:
>virtual_accelerator --debug  --sequences MEBT --bunch bunch_mebt_entrance_v0.dat --particle_number 10000 --refresh_rate 10

"""

import sys
import math
import time

from epics import pv as pv_channel
from epics import caget, caput

def calculateRMS(x_arr,y_arr):
	"""
	Calculates sigma and peak postion as x_avg for y = Function(x)
	"""
	x_avg = 0.
	y_sum = 0.
	for ind,x in enumerate(x_arr):
		y = y_arr[ind]
		x_avg += x*y
		y_sum += y
	x_avg /= y_sum
	x_sigma2 = 0.
	for ind,x in enumerate(x_arr):
		y = y_arr[ind]
		x_sigma2 += y*(x-x_avg)**2
	x_sigma2 /= y_sum
	x_sigma = math.sqrt(abs(x_sigma2))
	return (x_avg,x_sigma)

def wireScan(pos_start,pos_end,ws_speed,sleep_time = 1.0, ws_name = "04a"):
	"""
	Returns the sigma for horizontal and vertical WS distribution curves.
	pos_start,pos_end - positions in mm (like from -25. to 30.)
	ws_speed - defines the postition step 0.2, 0.3 or 0.5 mm
	sleep_time - time for VA respond - by default should be 1.0 sec.
	"""
	pos_set_pv = pv_channel.PV("MEBT_Diag:WS"+ws_name+":Position_Set")
	pos_pv = pv_channel.PV("MEBT_Diag:WS"+ws_name+":Position")
	speed_pv = pv_channel.PV("MEBT_Diag:WS"+ws_name+":Speed_Set")
	hor_signal_pv = pv_channel.PV("MEBT_Diag:WS"+ws_name+":Hor_Cont")
	ver_signal_pv = pv_channel.PV("MEBT_Diag:WS"+ws_name+":Ver_Cont")
	#----- Geometry of fork coeff
	coeff_fork = math.sqrt(2.0)/2
	#---- Retract fork -------------------
	speed_pv.put(100.)
	pos_set_pv.put(pos_start)
	time.sleep(2.0)
	#-------Set up speed and end point ---
	#-------This will start scan ---------
	speed_pv.put(ws_speed)
	pos_set_pv.put(pos_end)
	#-------------------------------------
	pos_arr = []
	hor_signal_arr = []
	ver_signal_arr = []
	pos = pos_pv.get()
	delta_pos = 0.01
	while(pos < pos_end - delta_pos):
		pos = pos_pv.get()
		hor_signal = hor_signal_pv.get()
		ver_signal = ver_signal_pv.get()
		time.sleep(sleep_time)
		print ("wire pos = %+6.3f"%pos," H,V signals = %12.5g %12.5g "%(hor_signal,ver_signal))
		pos_arr.append(pos)
		hor_signal_arr.append(hor_signal)
		ver_signal_arr.append(ver_signal)
	#---- Retract fork -------------------
	speed_pv.put(100.)
	pos_set_pv.put(pos_start)
	time.sleep(2.0)
	#-------------------------------------
	(x_avg,x_sigma) = calculateRMS(pos_arr,hor_signal_arr)
	(y_avg,y_sigma) = calculateRMS(pos_arr,ver_signal_arr)
	print ("debug (x_avg,x_sigma) = ",(x_avg,x_sigma*coeff_fork)," (y_avg,y_sigma) = ",(y_avg,y_sigma*coeff_fork))
	return (x_sigma*coeff_fork,y_sigma*coeff_fork)
	

#-------------------------------------------------------------------
#              START of the SCRIPT
#-------------------------------------------------------------------

sleep_time = 1.1

print ("================================")
#---- Let's set 1 RF Cavity in MEBT to zero amplitude
cav_amp_pv = pv_channel.PV("MEBT_LLRF:FCM1:CtlAmpSet")
print ("MEBT:FCM1 initial amplitude =",cav_amp_pv.get())
cav_amp_pv.put(0.)
time.sleep(sleep_time)
print ("MEBT:FCM1 new     amplitude =",cav_amp_pv.get())
print ("================================")
#---- Let's set quad #4 to zero gradient
quad_qv04_field_pv = pv_channel.PV(" MEBT_Mag:PS_QV04:B_Set")
print ("QV04 initial G[T/m] =",quad_qv04_field_pv.get())
quad_qv04_field_pv.put(0.)
time.sleep(sleep_time)
print ("QV04 new     G[T/m] =",quad_qv04_field_pv.get())
print ("================================")

quad_field_pv = pv_channel.PV(" MEBT_Mag:PS_QH03:B_Set")

quad_field_init = quad_field_pv.get()
print ("QH03 Initial G[T/m] =",quad_field_pv.get())

#---- Wire Scanner scan params. Positions in mm
pos_start = -25.0
pos_end = 30.0

#---- ws_speed defines step
ws_speed = 1.0
ws_sleep_time = 0.1

(hor_sigma,ver_sigma) = wireScan(pos_start,pos_end,ws_speed,ws_sleep_time)
print ("(hor_sigma,ver_sigma) =",(hor_sigma,ver_sigma))

#---- Let's define quad field scan
quad_field_start = 4.0 
quad_field_end   = 40.0
quad_field_step  = 4.0

#---------------------------
# To be continued
#---------------------------

sys.exit(0)

# Importing epics to use pyepics to talk to Epics
import epics

# Importing time to allow our script to sleep while we wait for changes to be registered.
import time

# Imports to help with plotting.
import numpy as np
import matplotlib.pyplot as plt

# Import curve fitting function to analyze Gaussian
from scipy.optimize import curve_fit

# PV names for the faraday cup we want to read.
fc_name = "FC"
# PV name for the BPM we want to read.
slit_name = "slit"

# Setup for the scan to make a plot later.
num = 50
slit_set_positions = np.linspace(21.5, 39, num)
slit_positions = np.zeros(num)
fc_charges = np.zeros(num)

# set slit speed
epics.caput(slit_name + ":Speed_Set", 3)
time.sleep(3)

print("Moving slit to desired initial position")
# move slit to desired initial position
epics.caput(slit_name + ":Position_Set", slit_set_positions[0])
while abs(epics.caget(slit_name + ":Position") - slit_set_positions[0]) > 1:
    time.sleep(1)

print("Start scan...")
for i in range(len(slit_set_positions)):
    # Set the slit position
    epics.caput(slit_name + ":Position_Set", slit_set_positions[i])

    # Sleep long enough for the slit to move.
    time.sleep(1)

    # Read the slit position and FC charge after update.
    slit_positions[i] = epics.caget(slit_name + ":Position")
    fc_charges[i] = epics.caget(fc_name + ":charge")

print("Scan finished.")

# Define the Gaussian function
def gauss(x, slw, x0, sigma):
    q0 = 10.1 # maximum charge with full beam
    return q0 * slw * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) / (np.sqrt(2*np.pi)*sigma)

print("Fitting beam profile to Gaussian distribution...")
# give an initial guess
guess = np.array([0.2, 30, 5])
parameters, covariance = curve_fit(gauss, slit_positions, fc_charges, guess)

fit_slw, fit_x0, fit_sigma = parameters

print("Fitting results:")
print("RMS beam size: " + str(fit_sigma) + " mm")
print("Center of beam: " + str(fit_x0) + " mm")
print("Slit width: " + str(fit_slw) + " mm")

fit_y = gauss(slit_positions, fit_slw, fit_x0, fit_sigma)
plt.plot(slit_positions, fc_charges, 'o', label='data')
plt.plot(slit_positions, fit_y, '--', label='fit')
plt.xlabel("Slit position [mm]")
plt.ylabel("FC charge [a.u.]")
plt.legend()
plt.show()

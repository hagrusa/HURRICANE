##########################################
    # File name: runProb1.py
    # Author: Harrison Agrusa, Ramsey Karim, Julian Marohnic
    # Date created: 11/12/2017
    # Date last modified: 11/19/2017
    # Description:  Create intial conditions files for problem 1, run code, produce plots
##########################################


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
#from astropy import constants as C
#import constants as cgs
#import time
import subprocess

####################################################################
############################## Set up ##############################
####################################################################


# No need for tree code -- 2 particles
theta_crit = 0.0
# No softening
epsilon = 0.0
# nsteps is on order of 5000 fyi
outfreq = 20
# Data is a massive dictionary that will be populated with all output from each run
data = {}
# Eccentricity
eccentricities = [0.5, 0.9]
# m = m1 = m2 = 1
m = 1
# Separation of particles
r = 1
velocities = [np.sqrt(2.*(1-e)) for e in eccentricities]
# 100 orbits
orbits = 100
# Orbital period
periods = [np.pi*np.sqrt(2./(1+e)**3) for e in eccentricities]
# Total integration times to get 100 orbits
times = [orbits*P for P in periods]
# Time steps
hs = [0.05, 0.003]
# Number of steps needed
Ns = [int(t/h) for t, h in zip(times, hs)]

# Change directories into the code directory, src/
os.chdir("src")
if len(sys.argv) == 2 and sys.argv[1] == "-m":
	subprocess.call("make")
elif len(sys.argv) > 1:
	print("You called the function with improper/too many arguments.")
	print("This function is normally called without arguments.")
	print("You may optionally add the flag '-m' to the end of the call to make/clean the C binaries.")
	sys.exit()

####################################################################
######################### Helper Functions #########################
####################################################################

def run_code(key):
	index, integrator = key.split('_')
	index = int(index[1]) - 1
	outdir_1 = "output1/" #output directory
	# Make output directory
	subprocess.call("mkdir {}".format(outdir_1), shell=True)
	subprocess.call("./nbodymain {} {} {} {} {} {} {} {}".format('prob1_{}.txt'.format(index),
		hs[index], Ns[index], epsilon, outfreq, integrator, theta_crit, outdir_1), shell=True)
	# For reference: nbody takes arguments of this form:
	# 	Usage: ./nbodymain data.txt h Nsteps epsilon outfreq integrator theta_critical output_dir starting_index(optional)
	outfiles = sorted(os.listdir(outdir_1))
	data[key] = np.empty((len(outfiles), 2, 7)) #2 particles, 7 data points
	# Get data
	for i in range(0, len(outfiles)):
		data[key][i, :, :] = np.genfromtxt(outdir_1 + str(outfiles[i]))
	# Remove directory, dont need it anymore
	subprocess.call("rm -r {}".format(outdir_1), shell=True)

def plot_xy(key):
	index, integrator = key.split('_')
	index = int(index[1]) - 1
	plt.plot(data[key][:, 0, 1], data[key][:, 0, 2], ".", label = r"$\rm Particle \ 1$")
	plt.plot(data[key][:, 1, 1], data[key][:, 1, 2], "x", label = r"$\rm Particle \ 2$")
	plt.xlabel(r"$\rm X \ Position$", fontsize = 18)
	plt.ylabel(r"$\rm Y \ Position$", fontsize = 18)
	plt.title(r"$\rm X \ vs \ Y \ e = {}$, {}".format(eccentricities[index], integrators[integrator]))
	plt.tight_layout()
	plt.legend(loc='best')
	plt.axis("equal")
	plt.savefig("XY_{}.pdf".format(key))
	plt.show()
	plt.close()

def plot_phase_energy(key):
	index, integrator = key.split('_')
	index = int(index[1]) - 1
	# r = sqrt(x^2+y^2+z^2)
	r_1 = np.sqrt(data[key][:, 0, 1]**2 + data[key][:, 0, 2]**2 + data[key][:, 0, 3]**2)
	Vr_1 = np.sum(data[key][:, 0, 4:7] * data[key][:, 0, 1:4])/r_1 #(v*r)/r

	r_2 = np.sqrt(data[key][:, 1, 1]**2 + data[key][:, 1, 2]**2 + data[key][:, 1, 3]**2)
	Vr_2 = np.sum(data[key][:, 1, 4:7] * data[key][:, 1, 1:4])/r_2

	plt.subplot(211)
	plt.plot(r_1, Vr_1, ".", label = r"$\rm Particle \ 1$")
	plt.plot(r_2, Vr_2, ".", label = r"$\rm Particle \ 2$")
	plt.xlabel(r"$\rm r \ (magnitude)$")
	plt.ylabel(r"$\rm Radial \ Velocity$")
	plt.title(r"$\rm Phase \ Diagram \ e = {}$, {}".format(eccentricities[index], integrators[integrator]))
	plt.legend()
	plt.tight_layout()
	plt.subplot(212)
	v1_sq = np.sum(data[key][:, 0, 4:7]**2, axis=1)
	v2_sq = np.sum(data[key][:, 1, 4:7]**2, axis=1)
	recip_r_sep = np.reciprocal(np.sqrt(np.sum((data[key][:, 1, 1:4] - data[key][:, 0, 1:4])**2, axis=1)))
	E = -recip_r_sep + (v1_sq + v2_sq)/2.
	E = (-2*E/(1 + eccentricities[index])) - 1
	time_array1 = np.arange(E.size) * hs[index]
	plt.plot(time_array1, E, '-', color='g')
	plt.title("Fractional Change in Total Energy vs Time")
	plt.ylabel("$\Delta E / E_{0}$")
	plt.xlabel("Time")
	plt.tight_layout()
	plt.savefig("phase_energy_{}.pdf".format(key))
	plt.show()
	plt.close()

####################################################################
################### Generate Initial Conditions ####################
####################################################################

def generate_input(index):
	f = open('prob1_{}.txt'.format(index), 'w')
	f.write("{} {} {} {} {} {} {}\n".format(m, -r/2.0, 0.0, 0.0, 0.0, -velocities[index]/2.0, 0.0))
	f.write("{} {} {} {} {} {} {}\n".format(m, r/2.0, 0.0, 0.0, 0.0, velocities[index]/2.0, 0.0))
	f.close()

for i in range(2):
	generate_input(i)

####################################################################
########################## Run Code ################################
####################################################################
sys.stdout.write("Running code.... \r")
sys.stdout.flush()

keys = ["e1_LF2", "e1_RK4", "e2_LF2", "e2_RK4"]
integrators = {"LF2": "Leapfrog", "RK4": "Runge-Kutta"}

for k in keys:
	run_code(k)

subprocess.call("rm prob1*.txt", shell=True)
if len(sys.argv) == 2 and sys.argv[1] == "-m":
	subprocess.call("make clean", shell=True)
os.chdir("..")


####################################################################
######################### Make Plots ###############################
####################################################################

for k in keys:
	plot_xy(k)
	plot_phase_energy(k)



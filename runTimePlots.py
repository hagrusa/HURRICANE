##########################################
    # File name: runTimePlots.py
    # Author: Harrison Agrusa, Ramsey Karim, Julian Marohnic
    # Date created: 11/12/2017
    # Date last modified: 11/20/2017
    # Description: create runtime vs number of particle plot
##########################################


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
#from astropy import constants as C
#import constants as cgs
import time
import subprocess

####################################################################
############################## Set up ##############################
####################################################################
"""
# Change directories into the code directory, src/
os.chdir("src")

if len(sys.argv) == 2 and sys.argv[1] == "-m":
    subprocess.call("make")
elif len(sys.argv) > 1:
    print("You called the function with improper/too many arguments.")
    print("This function is normally called without arguments.")
    print("You may optionally add the flag '-m' to the end of the call to make/clean the C binaries.")
    sys.exit()

os.chdir("..")
"""
def gen_particles(Nparticles):
    N = int(Nparticles)
    L = 100 #size of box --> V = (L)^3
    M = 1.0  #~earth size masses
    V = 0.0
    f = open('runtime_init.txt', 'w')
    for i in range(0, N):
        m = M
        x = (np.random.rand()-0.5)*L
        y = (np.random.rand()-0.5)*L
        z = (np.random.rand()-0.5)*L
        xdot = (np.random.rand()-0.5)*V
        ydot = (np.random.rand()-0.5)*V
        zdot = (np.random.rand()-0.5)*V
        if i != N-1: #dont want an empty last line
            f.write("{} {} {} {} {} {} {}\n".format(m, x, y, z, xdot, ydot, zdot))
        else:
            f.write("{} {} {} {} {} {} {}".format(m, x, y, z, xdot, ydot, zdot))
    f.close()







thetas = np.array([0.0, 0.2, 1.0])
#Nparticles = np.linspace(10, 1000, 10)
Nparticles = np.linspace(10, 1000, 10)

runtime = np.zeros((len(thetas), len(Nparticles)))
h = 0.01
Nsteps = 1000
epsilon = 1
outfreq = 1000000 #dont want to write output files
integrator = "LF2"



outdir = "output/" #output directory
init_cond = "runtime_init.txt"
subprocess.call("mkdir {}".format(outdir), shell=True)


for i in range(0, len(Nparticles)):
    gen_particles(Nparticles[i]) #generate initial conditions
    for j in range(0, len(thetas)):
        subprocess.call("mkdir {}".format(outdir), shell=True)
        start = time.time()
        subprocess.call("./nbodymain {} {} {} {} {} {} {} {}".format(init_cond,
            h, Nsteps, epsilon, outfreq, integrator, thetas[j], outdir), shell=True)
        stop = time.time()
        runtime[j, i] = stop - start
        subprocess.call("rm -r {}".format(outdir), shell=True)
    subprocess.call("rm {}".format(init_cond), shell = True)


plt.plot(Nparticles, runtime[0,:], 'o',label = "Theta = {}".format(thetas[0]))
plt.plot(Nparticles, runtime[1,:], 'o',label = "Theta = {}".format(thetas[1]))
plt.plot(Nparticles, runtime[2,:], 'o',label = "Theta = {}".format(thetas[2]))
plt.plot(Nparticles, (np.max(runtime[0,:])/ np.max(Nparticles**2)) * Nparticles**2, label = r"$N^{2}, \ (normalized)$")
plt.plot(Nparticles, (np.max(runtime[2,:])/ np.max(Nparticles*np.log10(Nparticles))) *  Nparticles*np.log10(Nparticles) , label = r"$NlogN, \ (normalized)$")
plt.xlabel(r"$\rm Number \ of \ Particles$", fontsize = 18)
plt.ylabel(r"$\rm Run \ Time \ [seconds]$", fontsize = 18)
plt.title(r"$\rm Tree \ Code \ Performance$", fontsize = 20)
plt.tight_layout()
plt.legend(loc='best')
plt.savefig("runtime.pdf")
plt.show()




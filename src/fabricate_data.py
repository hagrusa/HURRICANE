import numpy as np
import random as rand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as plt3d
import glob, os
import sys

inname = "data.txt"
outname = "treedata.txt"

def fabricate(num):
	file = open(inname, 'w')

	rand.seed(11)

	m = 128
	l = 0
	a, b = -m, m
	c, d = -l, l

	for i in range(num):
		line = "%f %f %f %f %f %f %f\n" % (rand.uniform(10000, 1000000),
			rand.uniform(a, b), 0, rand.uniform(a, b),
			0, rand.uniform(a, b), 0
			# rand.uniform(a, b), rand.uniform(c, d), rand.uniform(a, b),
			# rand.uniform(c, d), rand.uniform(a, b), rand.uniform(c, d)
		)
		file.write(line)
	file.close()

def check_tree():
	m, x, cx, y, cy, z, cz = np.loadtxt(inname, unpack = True)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	maxsize = 30
	plt3d.scatter(ax, x, y, z, c='r', s=maxsize*m/np.max(m), label="Particles")

	w, x, cx, y, cy, z, cz = np.loadtxt(outname, unpack = True)
	maxsize = 20
	# plt3d.scatter(ax, x, y, z, 'x', c='g', s=maxsize/4, label="Centers")
	# plt3d.scatter(ax, cz, cy, cz, c='hotpink', s=2*w/np.max(w), label="COMs")
	# for width, j, k, i in zip(w, x, y, z):
	for width, i, j, k in zip(w, x, y, z):
		# plt3d.scatter(ax, [i], [j], [k], c='r')
		for top in range(2):
			for side in range(2):
				plt3d.plot(ax, [i - width/2, i + width/2],
					[j - width/2 + side*width, j - width/2 + side*width],
					[k - width/2 + top*width, k - width/2 + top*width],
					c='k', linewidth=.2)
				plt3d.plot(ax, [i - width/2 + side*width, i - width/2 + side*width],
					[j - width/2, j + width/2],
					[k - width/2 + top*width, k - width/2 + top*width],
					c='k', linewidth=.2)
				plt3d.plot(ax, [i - width/2 + top*width, i - width/2 + top*width],
					[j - width/2 + side*width, j - width/2 + side*width],
					[k - width/2, k + width/2],
					c='k', linewidth=.2)
	plt.legend()
	plt.xlabel("X")
	plt.ylabel("Y")
	ax.set_zlabel("Z")
	plt.show()

def check_output():
	fname = "data.txt"
	m, x, vx, y, vy, z, vz = np.loadtxt(fname, unpack = True)
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	maxsize = 30
	plt3d.scatter(ax, x, y, z, c='k', s=maxsize*m/np.max(m), label="Before")
	plt.show()
	
def plot_results():
	os.chdir("result")
	all_data = []
	maxt = -1
	for f in sorted(os.listdir("./")):
		all_data.append(np.loadtxt(f, unpack=True))
		maxt += 1
	d = np.array(all_data)
	shape = d.shape
	m = d[0, 0, :]
	maxsize = 20
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	maxxyz = 0
	minsep = 1e8
	maxsep = 0
	for t in range(shape[0]):
		maxxyz = max(maxxyz, max(np.max(np.abs(d[t, 1, :])), np.max(np.abs(d[t, 2, :])), np.max(np.abs(d[t, 3, :]))))
	# 	if t == maxt:
	# 		print "at last t"
	# 		for particle1 in range(shape[2]):
	# 			for particle2 in range(shape[2]):
	# 				if particle1 == particle2:
	# 					continue
	# 				minsep = min(minsep, dist(d[t,:,particle1], d[t,:,particle2]))
	# 				maxsep = max(maxsep, dist(d[t,:,particle1], d[t,:,particle2]))
	# print minsep, maxsep
	for t in range(shape[0]):
		# plt.scatter(d[t, 1, :], d[t, 3, :], s=maxsize*m/np.max(m))
		plt3d.scatter(ax, d[t, 1, :], d[t, 2, :], d[t, 3, :], c='k', s=maxsize*m/np.max(m))
		plt.xlim([-maxxyz, maxxyz])
		plt.ylim([-maxxyz, maxxyz])
		ax.set_zlim([-maxxyz, maxxyz])
		fname = "graphic_result_%03d.png" % t
		sys.stdout.write("Writing "+fname+"\r")
		sys.stdout.flush()
		plt.savefig(fname)
		plt.cla()
		# plt.show()
	sys.stdout.write("\nDone!\n")
	sys.stdout.flush()

def dist(m_phase_array1, m_phase_array2):
	sep = m_phase_array2[1::2] - m_phase_array1[1::2]
	sep = sep*sep
	sep = np.sqrt(np.sum(sep))
	return sep


if len(sys.argv) > 1:
	if sys.argv[1] == "fabricate":
		fabricate(int(sys.argv[2]))
	elif sys.argv[1] == "check":
		check_output()
else:
	plot_results()

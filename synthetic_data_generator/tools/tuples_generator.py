#!/usr/bin/python
# 
# Program: Script to generate a dataset of points ordered by their application timestamps.
#		   The offset between the application timestamps of two consecutive tuples (inter-arrival
#          time) is generated using an exponential distribution if the dispersion index is 0,
#          otherwise we use a MAP(2) for generating bursty arrivals with a given index of dispersion.
# 
#          VERSION REVISED BASED ON NUMPY ARRAYS (FASTER THAN THE ORIGINAL SCRIPT)
# Date:	   November 2016
# Author:  Gabriele Mencagli <mencagli@di.unipi.it>

import sys
import getopt
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from burst import MAP_Distribution

# top level names (configuration parameters)
corr = 0. # correlation factor, corr must be in the range [-1,1]
val_stddiv = 30. # stddev of the points' coordinates
val_mean = 500. # mean of the points' coordinates

# print the arrivals per second into a text file
def print_arrivals(arrivals):
	file_out = open("arrivals.txt", 'w')
	count = 0
	for i in range (0, len(arrivals)):
		file_out.write(str(count) + ' ' + str(arrivals[i]) + '\n')
		count = count + 1
	file_out.close()

# plot the arrivals per second
def plot_arrivals(arrivals):
	ids = range (0, len(arrivals))
	plt.plot(ids, arrivals)
	plt.show()

# function to return the element at position i,j of the correlation matrix
def fij(i, j):
	if(i!=j): return corr
	else: return 1.

# function to check if a matrix is positive semi-definite
def isPSD(A, tol=1e-8):
	E,V = linalg.eigh(A)
	return np.all(E > -tol)

# function to adjust a correlation matrix in order to be likely positive semi-definite
def correctMatrix(M):
	D,V = linalg.eig(M)
	V1 = V[:,1]
	M2 = M + V1*V1.T*(-1*np.spacing(np.single(D[1].real))-D[1])
	return M2.real

# main function
def main():
	global corr # in this function corr can be modified
	if(len(sys.argv) !=  13):
		sys.stdout.write("Usage: " + sys.argv[0] + " -n <no. tuples> -d <no. dimensions> -c <correlation [-1,1]> -r <no. tuples per second> -i <dispersion index> -f <output filename>" + '\n')
		sys.exit(1)
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'n:d:c:r:i:f:')
	except getopt.GetoptError:
		sys.stdout.write("Usage: " + sys.argv[0] + " -n <no. tuples> -d <no. dimensions> -c <correlation [-1,1]> -r <no. tuples per second> -i <dispersion index> -f <output filename>" + '\n')
		sys.exit(1)
	for opt, arg in opts:
		if opt in ('-n'):
			n = int(arg) # number of tuples to generate
		elif opt in ('-d'):
			d = int(arg) # number of dimensions per tuple
		elif opt in ('-c'): # correlation factor of tuple's attributes
			if(float(arg) > -1 and float(arg) < 1):
				corr = float(arg)
			else:
				sys.stdout.write("Invalid correlation factor: using c=0" + '\n')
				corr = 0.
		elif opt in ('-r'):
			rate = float(arg) # rate per second in no. of tuples
		elif opt in ('-i'):
			disp_index = float(arg) # dispersion index (0 no bursts, >0 bursty distribution)
		elif opt in ('-f'):
			output_file = str(arg) # output file of the generated dataset
		else:
			sys.stdout.write("Usage: " + sys.argv[0] + " -n <no. tuples> -d <no. dimensions> -c <correlation [-1,1]> -r <no. tuples per second> -i <dispersion index> -f <output filename>" + '\n')
			sys.exit(1)
	# different point distributions
	if(corr != 0): # anticorrelated and correlated cases
		means = map(lambda x: val_mean, xrange(d))
		# compute the correlation matrix
		corrM = np.fromfunction(np.vectorize(fij), (d, d), dtype=float)
		# check whether the correlation matrix is positive semi-definite or not
		if(not isPSD(corrM)):
			sys.stdout.write("Correlation matrix is not positive semi-definite, I try to correct it..." + '\n')
			corrM = correctMatrix(corrM)
			if(not isPSD(corrM)):
				sys.stdout.write("Correction failed!" + '\n')
				sys.exit(1)
			else: sys.stdout.write("Correction done!" + '\n')
		# compute the covariance matrix
		covs = corrM * (val_stddiv**2)
		m = np.random.multivariate_normal(means, covs, n)
	else: # independent case
		m = (np.random.rand(n, d) * (val_mean*2))
	# create the dataset as a numpy array
	dataset = np.zeros((n, d+1))
	dataset[:,:-1] = m
	# average inter-arrival time in timeunits (usec)
	TA = 1000000.0/rate
	# if disp_index>0 we use the MAP distribution to generate the inter-arrival times
	if(disp_index > 0):
		distribution = MAP_Distribution(TA, disp_index);
	# variables for the application timestamps generation
	base = 0
	curr_incr = 0
	avg_ta = 0.
	id_tuple = 0.
	# variables for computing the arrivals per second
	cnt_tuples = 0 # no. of tuples transmitted in the current second
	threshold = 1000000 # threshold to detect the end of the next second
	arrivals_sec = [] # counts of arrivals per second
	# for all the tuples in the dataset
	for i in range(n):
		dataset[i][d] = base; # application timestamp of the i-th tuple
		# update the average of the inter-arrival time
		id_tuple = id_tuple + 1
		avg_ta = avg_ta + (1/id_tuple) * (curr_incr - avg_ta)
		# update the measurements for computing the arrivals per second
		if(base > threshold):
			threshold = threshold + 1000000
			arrivals_sec.append(cnt_tuples)
			cnt_tuples = 1
		else:
			cnt_tuples = cnt_tuples + 1
		# if disp_index>0 we use a MAP(2) distribution for the inter-arrival times
		if(disp_index > 0):
			curr_incr = distribution.next_sample()
		# otherwise we use an exponential distribution
		else:
			curr_incr = np.random.exponential(TA)
		base += curr_incr
	# save the dataset in a textual file
	np.savetxt(output_file, dataset, fmt='%1.9f')
	sys.stdout.write("Mean inter-arrival time: " + str(avg_ta) + ' usec\n')
	sys.stdout.write("Mean arrival rate: " + str(1000000/avg_ta) + ' tuples/sec\n')
	print_arrivals(arrivals_sec)
	#plot_arrivals(arrivals_sec)

if __name__ ==  "__main__":
    main()

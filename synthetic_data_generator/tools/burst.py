#!/usr/bin/python
# 
# Program: Script to generate a MAP(2) distribution to obtain bursty
#		   application timestamps in a dataset.
# Date:    March 2016
# Author:  Gabriele Mencagli <mencagli@di.unipi.it>

import sys
import random

# class definition of a MAP(2) distribution
class MAP_Distribution:
	mean = 0. # mean of the distribution
	d_index = 0. # index of dispersion
	mean_fast = 0. # mean in the fast state
	mean_slow = 0. # mean in the slow state
	pfs = 0. # prob transition fast state -> slow state
	psf = 0. # prob transition slow state -> fast state
	active_state = 'slow' # initially we start in the slow state

	# constructor
	def __init__(self, _mean, _d_index):
		self.mean = _mean
		self.d_index = _d_index
		self.fitting() # start the fitting

	# method to obtain a new sample
	def next_sample(self):
		# calculate the sample
		sample = 0.
		if(self.active_state == 'slow'):
			# slow state case
			sample = random.expovariate(1/self.mean_slow)
		elif(self.active_state == 'fast'):
			# fast state case
			sample = random.expovariate(1/self.mean_fast)
		# update the state of the MAP(2)
		r =  random.uniform(0, 1) # random number in [0.1)
		if((self.active_state == 'slow') and (r <= self.psf)):
			self.active_state = 'fast'
		elif((self.active_state == 'fast') and (r <= self.pfs)):
			self.active_state = 'slow'
		return sample

	# method to compute the index of dispersion
	def get_dindex(self):
		rate_fast = 1/self.mean_fast
		rate_slow = 1/self.mean_slow
		numerator = 2*self.pfs*self.psf*((rate_fast-rate_slow)**2)
		denominator = (self.pfs+self.psf)*((rate_fast*self.pfs+rate_slow*self.psf)**2)
		return 1+(numerator/denominator)

	# method to compute the lag-1 autocorrelation
	def get_lag1(self):
		val = 2*((self.psf*self.mean_fast**2)/(self.psf+self.pfs)) + 2*((self.pfs*self.mean_slow**2)/(self.psf+self.pfs))
		lag = 0.5*(1-self.psf-self.pfs)*(1-((self.mean**2)/(val-self.mean**2)))
		return lag

	# fitting procedure for the MAP(2) distribution
	def fitting(self):
		F = 1.0001 # scaling parameter
		ACF = 0. # lag-1 autocorrelation
		threshold = 0.4 # threshold for lag-1 autocorrelation
		sys.stdout.write('Starting fitting...' + '\n')
		# while loop until the lag-1 becomes greater than the threshold
		while ACF < threshold:
			self.mean_slow = F*self.mean
			self.mean_fast = self.mean/F
			self.pfs = 0.01
			self.psf = self.pfs*((self.mean_slow-self.mean)/(self.mean-self.mean_fast))
			MAXITER = 1000; # maximum number of iterations
			# while loop until the index of dispersion is within 1% of the target one
			while(MAXITER > 0 and abs(self.get_dindex()-self.d_index)>0.01*self.d_index):
				MAXITER = MAXITER-1
				if(self.get_dindex() > self.d_index):
					# the index of dispersion is high so increase pfs to be in the slow mode longer
					self.pfs = self.pfs/random.uniform(0, 1)
					self.psf = self.pfs*((self.mean_slow-self.mean)/(self.mean-self.mean_fast))
				elif(self.get_dindex() < self.d_index):
					self.pfs = self.pfs*random.uniform(0, 1)
					self.psf = self.pfs*((self.mean_slow-self.mean)/(self.mean-self.mean_fast))
			self.d_index = self.get_dindex()
			ACF = self.get_lag1()
			F = F*2;
		sys.stdout.write('...complete' + '\n')
		sys.stdout.write('Index of dispersion of the MAP: ' + str(self.get_dindex()) + '\n')

if __name__ ==  "__main__":
	if(len(sys.argv) !=  4):
		sys.stdout.write("Usage: " + sys.argv[0] + " <mean> <dispersion index> <no. points>" + '\n')
		sys.exit(0)
	mean = float(sys.argv[1])
	index = float(sys.argv[2])
	points = int(sys.argv[3])
	# example of usage of the distribution
	dist = MAP_Distribution(mean, index)
	avg = 0.
	count = 0.
	for i in range(0, points):
		count = count + 1
		sample = dist.next_sample()
		sys.stdout.write(str(sample) + '\n')
		avg = avg + (1/count) * (sample - avg)
	sys.stdout.write("Average: " + str(avg) + '\n')

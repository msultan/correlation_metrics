#!/usr/bin/env python
from __future__ import division
import posMutualCode
from posMutualCode import * 
import numpy as np
def testCodeWrapper():
	import matplotlib.pyplot as mlp


	tN=50
	N=800
	k=5
	mean=[3,2]
	cov=[[1,0.9],[0.9,1]]
	data=np.zeros((N,2))
	
	tD=np.zeros(tN)
	

	for t in range(tN):
		#get random data	
		data=np.random.multivariate_normal(mean,cov,N)
		mutual=mutual_nearest_neighbors(N,k,data)
		tD[t]=mutual
	 
	det=np.linalg.det(cov)
	print "Analytical Solution:",-0.5*np.log(det)
	print 'Average of %d Trial for %d datapoints is %f with standard deviation %f'%(tN,N, np.mean(tD),np.std(tD))
	mlp.plot(tD)
	mlp.axhline(-0.5*np.log(det))
	mlp.show()

	return

def main(assignFile,projectFile,test):
	'''
	Mutual information Calculator for the positional Vectors of a specified residues. This code is based
	of the work of Kraskov,McClendon and Lange. 
	Parameters:
	----------
	assignment File: File with the macro Assignments
	project File: The project file 
	
	iterations: how many iterations/permutations for each data
	alignment File:

	align_indices : np.ndarray or None 
						atom indices to use in the alignment step
	atom_indices : np.ndarray or None
					atom indices to use when calculating distances
	Output:
	----------
	multiple *.dat files which has mutual information for each state in the assignments file
	'''
	import msmbuilder as m
	import numpy as np
	#load the assignments
	macroAssignments=m.io.loadh(assignFile)
	#get the actual assignment
	macroAssignments=macroAssignments['arr_0']

	macroAssignmentsMax=np.max(macroAssignments)

	for i in range(0,macroAssignmentsMax+1):
		print "Currently Calculating Mutual for Start %d",i
		#macroAssignmentsIndex=np.where(macroAssignments==i)


	return


def parse_commandline():
	import os
	import optparse
	parser = optparse.OptionParser()
	parser.add_option('-a', '--assignmentFile', dest='assignFile', help='Assignment File')
	parser.add_option('-p','--projectFile',dest='projectFile',help='Project File')
	parser.add_option('-i', '--total_iterations', dest='i', type="int", default=1, help=' total iterations')
	parser.add_option('--test', dest='test', default=0, help='test the code')
	parser.add_option('-s', '--skip_rows', dest='s',type='int', default=0,help='how many rows to skip')
	parser.add_option('-b', '--bins',dest='bin_n', type='int', default=24, help='The number of Bins used to bin data(Default 24). Too few bins can lead to problems')
	parser.add_option('-w', '--windows',dest='numWin', type='int', default=1, help='Whether or not to use windows in the calculation')
	(options, args) = parser.parse_args()
	return (options, args)

if __name__ == '__main__':
	(options, args) = parse_commandline()
	if (options.test): 
		testCodeWrapper()
	else:
		main(assignFile=options.assignFile,projectFile=options.projectFile)
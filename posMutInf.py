#!/usr/bin/env python

def mutualWrapper:
	try:
		import posMutualCode
	except ImportError:
		raise Exception("posMutualCode NOT FOUND ")

	N=len(data(:,1))
	mutual_nearest_neighbors(N,k,data)

def main(assignFile):
	import msmbuilder as m
	import numpy as np
	#load the assignments
	macroAssignments=m.io.loadh(assignFile)
	#get the actual assignment
	macroAssignments=macroAssignments['arr_0']

	macroAssignmentsMax=np.max(macroAssignments)

	for i in range(0,macroAssignmentsMax+1):
		print "Currently Calculating Mutual for Start %d",i
		macroAssignmentsIndex=np.where(macroAssignments==i)


	return


def parse_commandline():
    import os
    parser = optparse.OptionParser()
    parser.add_option('-aF', '--assignmentFile', dest='assignFile', type="file", help='Assignment File')
    parser.add_option('-pF','--projectFile',dest='projectFile',type="File",help='Project File')
    parser.add_option('-i', '--total_iterations', dest='i', type="int", default=1, help=' total iterations')
    parser.add_option('--test', dest='test', default=0, help='test the code')
    parser.add_option('-s', '--skip_rows', dest='s',type='int', default=0,help='how many rows to skip')
    parser.add_option('-b', '--bins',dest='bin_n', type='int', default=24, help='The number of Bins used to bin data(Default 24). Too few bins can lead to problems')
    parser.add_option('-w', '--windows',dest='numWin', type='int', default=1, help='Whether or not to use windows in the calculation')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == '__main__':
	(options, args) = parse_commandline()
	main(assignFile=options.assignFile,projectFile=options.projectFile)
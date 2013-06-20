#!/usr/bin/env python
from __future__ import division
import posMutualCode
from posMutualCode import * 
from testMutualCode import testPositionMutualInformationCode
import numpy as np

def positionalMutualCalculator(assignFile,projectFile,gensFile,atomIndices,states):
	'''
	Mutual information Calculator for the positional Vectors of a specified \
	residues. This code is based
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
	multiple *.dat files which has mutual information for \
	each state in the assignments file
	'''
	import msmbuilder as m
	from msmbuilder import Trajectory
	import numpy as np
	import lprmsd
	import os
	from collections import defaultdict
	from IPython import parallel
	#setting up the MAP
	client_list=parallel.Client(profile='mpi')
   	print "Running on:",len(client_list.ids)
  	view = client_list.load_balanced_view()

	#Load the Atom Indices 
	atomIndices=np.loadtxt(options.atomIndices,np.int)
	#making a dictionary for fast access to location of where the final value 
	#will end up in the matrix
	atomDict={}
	for i in atomIndices:
		atomDict[atomIndices[i]]=i
	# Load the project 
	prj = m.Project.load_from(projectFile)
	#load the assignments
	macroAssignments = m.io.loadh(assignFile)
	#get the actual assignment
	macroAssignments = macroAssignments['arr_0']

	macroAssignmentsMax = np.max(macroAssignments)


	#eventually Need to update this so that only certain states are tabulated

	if -1 == states:
		print "Calculating Mutual Information for all states"
		#currently going to calculate MI for all states
		states = np.arange(macroAssignmentsMax+1)

	#setting up Lee Ping's Metric which None for permute indices and \
	#alternative atom indices
	rmsd_metric=lprmsd.LPRMSD(atomIndices,None,None)

	#loading the generator file and creating a trajectory out of it.
	if os.path.exists(gensFile):
		gensTraj = Trajectory.load_trajectory_file(gensFile)
		pgenTraj = rmsd_metric.prepare_trajectory(gensTraj)



	#creating an inverse assignment dictionary to save all \
	#frames from all trajectories to a single 
	stateAssignmentDict=defaultdict(lambda:[])

	#{key:value} where key is the state and value is a \
	#list of tuple where each tuple has form(trjIndex,frmIndex)

	for trjIndex in xrange(macroAssignments.shape[0]):
		for frmIndex in xrange(macroAssignments.shape[1]):
			stateAssignmentDict[macroAssignments[trjIndex,frmIndex]]\
			.append((trjIndex,frmIndex))

	#number of neighbor
	k=5

	#loop through the states
	for s in states:
		mMat=np.zeros((len(atomIndices),len(atomIndices)))
		print "Calculating MI for state %s"%s
		if len(stateAssignmentDict[s])==0:
			raise ValueError('No Assignments to state %s'%s)
		#getting all conformation
		confs=stateAssignmentDict[s]
		
		#creating a frame dictionary so that i can pull those.
		FrameDict = {}
		for (traj, frame) in confs:
			FrameDict.setdefault(traj,[]).append(frame)
		
		#getting an empty traj
		clusterTraj=prj.empty_traj()
		#getting a set of what trajectories we need to query
		TrajNums=set(i[0] for i in confs)

		#getting only the frames we want for this state
		for currTrj in TrajNums:
			T=prj.load_traj(currTrj)[np.array(FrameDict[currTrj])]
			clusterTraj += T

		print "Loaded %i conformations"%len(clusterTraj)
		#Now, we should have clusterTraj, we can prepare it
		pclusterTraj=rmsd_metric.prepare_trajectory(clusterTraj)
		rmsd,xout=rmsd_metric.one_to_all_aligned(pgenTraj, pclusterTraj, s)
		#xout is the aligned trajectory, we need to subtract every value in it 
		#from the generator to the deviation from the mean
		N=len(xout)
		randomT=np.random.randint(N)
		randomI=np.random.randint(len(xout[0]))
		sanityTest=xout[randomT,randomI]
		
		#doing the actual subtraction
		xout=xout-pgenTraj[s]['XYZList']
		


		#simple test, basically subtract the ensemble average from a random 
		#atom index at a random time step and see if they are equal. 
		sanityTestValue=(sanityTest-pgenTraj[s]['XYZList'][0,randomI])
		assert(((xout[randomT,randomI]) == (sanityTestValue)).all())

		jobs=[]
		#for the positional vectors
		for indexTracker,atomindexI in enumerate(atomIndices):
			for indexTracker2,atomindexJ in enumerate(atomIndices[indexTracker:]):
				data=np.zeros((N,6))
				data[:,:3]=xout[:,atomindexI]
				data[:,3:]=xout[:,atomindexJ]
				job=(N,k,atomindexI,atomindexJ,data)
				jobs.append(job)
				#mMat[indexTracker][indexTracker2] = mutual_nearest_neighbors(N,k,data)

		results=view.map(mutual_nearest_neighbors,*zip(*jobs))
		all_mutuals = results.get()
        	for i,job in enumerate(jobs):
        		print atomDict[results[i][0]],results[i][1]
        		mMat[atomDict[results[i][0]]][atomDict[results[i][1]]]=\
        		mMat[atomDict[results[i][1]]][atomDict[results[i][0]]]=\
        		results[i][-1]

		np.savetxt('%s.dat'%s,mMat)

	return 0


def parse_commandline():
	import os
	import optparse
	parser = optparse.OptionParser()
	parser.add_option('--test', dest='test', default=0, \
					help='test the code')
	parser.add_option('-a', '--assignmentFile',default='Macro4/MacroAssignments.h5'\
					,dest='assignFile',help='Assignment File')
	parser.add_option('-p','--projectFile',dest='projectFile',\
					help='Project File',default='ProjectInfo.h5')
	parser.add_option('-g','--genFile',default='Data/Gens.lh5',\
					help='The Generator File that defines the "ensemble average"')
	parser.add_option('-i', '--total_iterations', dest='i', type="int", \
					default=1, help=' total iterations')
	parser.add_option('-s', '--states',dest='states', nargs="+", \
					default=-1, help=' States to Use to Calculate MI')
	parser.add_option('--atomIndices',  default='AtomIndices.dat',\
					help='atomIndices to align against and calculate the MI for')
	parser.add_option('-w', '--windows',dest='numWin', type='int', default=1,\
					help='Whether or not to use windows in the calculation')
	(options, args) = parser.parse_args()
	return (options, args)

if __name__ == '__main__':
	(options, args) = parse_commandline()


	if (options.test): 
		testPositionMutualInformationCode()
	else:
		positionalMutualCalculator(options.assignFile,options.projectFile,options.genFile,\
			options.atomIndices,options.states)




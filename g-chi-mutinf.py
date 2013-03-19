#!/usr/bin/env python
from __future__ import division
import numpy  # for histogram stuff
import optparse
import os
import glob
import itertools
import random
import multiprocessing
import pickle  # for compressing the data not really needed right now
from IPython import parallel
from IPython.parallel import require
from IPython.parallel import Client
from collections import defaultdict
import sys
import scipy
import string
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> c88e3732e09277dccd44013b3520f720794f6a57




import time                                                
>>>>>>> 73f072588370ee8f55c43f69e4a08ea4edd47e2f

def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print 'Ran %r residue pairs in %2.2f sec' % \
              (options.t, te-ts)
        return result

    return timed


def mutual_information_from_files(res_name1, res_id1, res_name2, res_id2,			shuffle, dir, skiprows,bin_n, test):

	#mutual_information_from_files
	"""
	Contact: msultan at stanford dot edu

	The g-chi-mutinf.py utility is designed to calculate the mutual information between the dihedral angles of protein amino acid in an attempt to understand it mechanistic pathway
	Args:
		dir: the directory with the dihedral angles(with extension .xvg). This program is designed to work with the g-chi gromacs utility Please Make sure that the dihedral files are sequentially labelled(1-> total_n_residues)
		total_n_residues: The total number of residues. 
		n_iterations: How many items the data needs to be scrambled
		skiprows: The number of rows that need to be skipped in each .xvg file(eg g-chi out put has the first 12 rows as junk data) 
		bin_n: number of bins
		test: Boolean for testing the script
	Returns:
		Mutual Information
	"""
	import numpy  #for histogram stuff
	import os
	import glob
	from scipy.stats.stats import pearsonr
	if not test:
		res_id1_range=[-180,180]
		res_id2_range=[-180,180]
		# example file namephiALA1.xvg.all_data_450micro

		filename_id1 = glob.glob('%s%s???%d.*' %(dir,res_name1,res_id1))
		filename_id2 = glob.glob('%s%s???%d.*' %(dir,res_name2,res_id2))
		if	filename_id1 and filename_id2:
			dihedral_id1 = numpy.loadtxt(filename_id1[0],usecols=(1,),skiprows=skiprows)	  # getting just the angles from the filepointer file
			dihedral_id2 = numpy.loadtxt(filename_id2[0],usecols=(1,),skiprows=skiprows)	 # getting just the angles from the filepointer file
			if shuffle:
				dihedral_id1 = numpy.random.permutation(dihedral_id1)
				dihedral_id2 = numpy.random.permutation(dihedral_id2)
		else:
			return 0.0
	if test:
		res_id1_range=[-3,3]
		res_id2_range=[-3,3]
		dihedral_id1 = numpy.random.normal(loc=0.0,scale=1.0,size=10000)
		dihedral_id2 = numpy.random.normal(loc=0.0,scale=1.0,size=10000)
#The magic begins here:The mutual information is defined as the 2d	sum of joint probability of x,y multiplied by the log(joint/marginal probabilites probabilities)

	
	(joint_pdf_bin_count,id1_bin_edges,id2_bin_edges) = numpy.histogram2d(dihedral_id1,dihedral_id2,bins=bin_n,range=[res_id1_range, res_id2_range],normed=False)
	
	dihedral_id1_count=numpy.sum(joint_pdf_bin_count,axis=0)
	dihedral_id2_count=numpy.sum(joint_pdf_bin_count,axis=1)
	
	dx=id1_bin_edges[1]-id1_bin_edges[0]
	N=len(dihedral_id1)
	H_x=-numpy.nansum((dihedral_id1_count/N) *numpy.log(dihedral_id1_count/N))+numpy.log(dx)
	H_y=-numpy.nansum((dihedral_id2_count/N) *numpy.log(dihedral_id2_count/N))+numpy.log(dx)
	H_x_y=-numpy.nansum((joint_pdf_bin_count/N) *numpy.log(joint_pdf_bin_count/N))+numpy.log(dx*dx)
	
	mutual=H_x+H_y-H_x_y
	
	if not test:
		return mutual
	if test:
		print H_x
		det=numpy.linalg.det(numpy.cov(dihedral_id1,dihedral_id2))
		print "The analytical solution for joint entropy is ",(0.5*numpy.log(det * (2*numpy.pi*numpy.e)**2))
		print "The calculated joint entropy is",H_x_y
		print "The Mutual Information should be", ((0.5*numpy.log(2*numpy.pi*numpy.e))+(0.5*numpy.log(2*numpy.pi*numpy.e)) -(0.5*numpy.log(det * (2*numpy.pi*numpy.e)**2))), "and the number calculated is",mutual
#		
@timeit
def main(dir,total_n_residues,n_iterations,skiprows,bin_n, test):
	olderr = numpy.seterr(all='ignore') 
	c=Client(profile='default')
	dihedral_names = ['chi1','chi2','chi3','phi','psi']
	jobs = []
	if test:
		print "TESTING BEGINNING:"
		print"Testing the System by calculating the entropy and mutual information for two normal distributions centred around 0 and 1 with std. dev. of1.The entropy for a single distribution should analyticaly be", (0.5*numpy.log(2*numpy.pi*numpy.e)),"it is calculated as"
		mutual_information_from_files('junk',1,'junk',1,True,dir,skiprows,bin_n,True)
		print "TESTING FINISHED"
		sys.exit()
	else:
		for id1 in range(1, total_n_residues+1):
			for id2 in range(id1, total_n_residues+1):
				for ic in range(0, n_iterations):
					for i, name1 in enumerate(dihedral_names):
						for name2 in dihedral_names[i:]:
							# construct the jobs
							if ic == 0:
								job = (name1, id1, name2, id2, False,dir,skiprows,bin_n,False)
							else:
								job = (name1, id1, name2, id2, True,dir,skiprows,bin_n,False)
							jobs.append(job)
	print "Running on:",len(c.ids)
	view = c.load_balanced_view()
                #async_results = []
                #for i, job in enumerate(jobs):
# ar = view.apply_async(mutual_information_from_files,*job)
                # async_results.append(ar)
                #print "Submitted:",len(async_results),"Jobs"
                #c.wait(async_results)
                #all_mutuals=[ar.get() for ar in async_results]
	result = view.map_async(mutual_information_from_files, *zip(*jobs))
	result.wait()
	all_mutuals = result.get()
	grids = {}
		
		final_grid = numpy.zeros(((total_n_residues),(total_n_residues)))
		average_grid = numpy.zeros(((total_n_residues),(total_n_residues)))
		for i,job in enumerate(jobs):
			name1 = job[0]
			id1 = job[1]
			name2 = job[2]
			id2 = job[3]
			shuffle = job[4]
			if shuffle:
				average_grid[id1-1][id2-1] += all_mutuals[i]
				average_grid[id2-1][id1-1] += all_mutuals[i]
			else:
				final_grid[id1-1][id2-1] += all_mutuals[i]
				final_grid[id2-1][id1-1] += all_mutuals[i]
		
		final_grid=final_grid-(average_grid/n_iterations)		
# 		for i, job in enumerate(jobs):
# 			name1 = job[0]
# 			id1 = job[1]
# 			name2 = job[2]
# 			id2 = job[3]
# 			shuffle = job[4]
# 			if shuffle:
# 				grid_name = 'grid_ave_'+ (name1) + '_' +(name2)
# 				if grid_name not in grids:
# 					grids[grid_name] = numpy.zeros(total_n_residues*total_n_residues)
# 				grids[grid_name][((id1-1)*total_n_residues)+id2-1] += all_mutuals[i]
# 			else:
# 				grid_name = 'grid_'+ (name1) + '_' +(name2)
# 				if grid_name not in grids:
# 					grids[grid_name] = numpy.zeros(total_n_residues*total_n_residues)
# 				grids[grid_name][((id1-1)*total_n_residues)+id2-1] = all_mutuals[i]
# 
# #Notes
# #the [((id1-1)*total_n_residues)+id2-1] index compresses a 2d array into a 1d array
# #for example, 1,1 is stored in (0*n+0)
# #1,n is stored in (0+n-1)
# #n,n is stored in(n*(n-1)+n-1)
# Random text twice now
# 
# 		final_grid=numpy.zeros(((total_n_residues),(total_n_residues)))
# 		for i, name1 in enumerate(dihedral_names):
# 			for name2 in dihedral_names[i:]:
# 				grid_name_1 = 'grid_'+ name1+ '_' +name2
# 				grid_name_2 = 'grid_ave_'+ name1 + '_' +name2
# 				for id1 in range(1, total_n_residues):
# 					for id2 in range(id1, total_n_residues):
# 						if n_iterations > 1:
# 							final_grid[id1-1][id2-1] += grids[grid_name_1][((id1-1)*total_n_residues)+id2-1] - grids[grid_name_2][((id1-1)*total_n_residues)+id2-1]/n_iterations
# 							final_grid[id2-1][id1-1] = final_grid[id1-1][id2-1]
# 						else:
# 							final_grid[id1-1][id2-1] += grids[grid_name_1][((id1-1)*total_n_residues)+id2-1]
# 							final_grid[id2-1][id1-1] = final_grid[id1-1][id2-1]
		file=open('%s'%dir+'mut-inf-mat.txt', 'w')
		numpy.savetxt(file, final_grid)


def parse_commandline():
	parser = optparse.OptionParser()
	parser.add_option('-d', '--directory', dest='dir',default="../", help='directory with chi1,phi,psi files')
	parser.add_option('-t','--total_residues',dest='t',type="int",help='total		number of residues')
	parser.add_option('-i', '--total_iterations', dest='i', type="int", default=1, help='	total iterations')
	parser.add_option('--test', dest='test', default=0, help='test the code')
	parser.add_option('-s', '--skip_rows', dest='s',type='int', default=0,help='how many rows to skip')
	parser.add_option('-b', '--bins',dest='bin_n', type='int', default=24, help='The number of Bins used to bin data(Default 24). Too few bins can lead to problems')
	(options, args) = parser.parse_args()
	return (options, args)

#function main if namespace is main

if __name__ == "__main__":
	(options, args) = parse_commandline()
	main(dir=options.dir, total_n_residues=options.t,n_iterations=options.i,skiprows=options.s,bin_n=options.bin_n, test=options.test)

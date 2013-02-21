#!/bin/python
from __future__ import division
import numpy #for histogram stuff
import optparse
import pylab
import os
import glob
import itertools
import pprint
import random
import multiprocessing
import pickle #for compressing the data not really needed right now
from IPython import parallel
from IPython.parallel import require
from IPython.parallel import Client
from collections import defaultdict

client = parallel.Client(profile='default') 
lbview = client.load_balanced_view()
dview = client[:]


def mutual_information_from_files(res_name1, res_id1, res_name2, res_id2, shuffle,dir):
	import numpy #for histogram stuff
	import optparse
	import pylab
	import os
	import glob
	import itertools
	import pprint
	import random
	import multiprocessing

	#return 'res_name1=%s res_id1=%s res_name2=%s res_id2=%s, shuffle=%s' % (res_name1, res_id1, res_name2, res_id2, shuffle)
	string_residue_counter_i = str(res_id1)
	string_residue_counter_j = str(res_id2)
	filename_i = glob.glob('%s' %dir +'%s' %res_name1 + '???' + string_residue_counter_i+ '.xvg.all_data_450micro')
	filename_j = glob.glob('%s' %dir +'%s' %res_name2 + '???' + string_residue_counter_j+ '.xvg.all_data_450micro')
	if  filename_i and filename_j:
		filepointer_i = open(filename_i[0],"r")
		phi_i = numpy.loadtxt(filepointer_i,usecols=(1,)) # getting just the angles from the filepointer file 
		filepointer_j= open(filename_j[0],"r")
		phi_j = numpy.loadtxt(filepointer_j,usecols=(1,)) # getting just the angles from the filepointer file
		if shuffle:
			phi_i=numpy.random.permutation(phi_i)
			phi_j=numpy.random.permutation(phi_j)
#The magic begins here 
#the mutual information is defined as the 2d sum of joint probability of x,y multiplied by the log(joint/individual probabilities)
		(H,x,y)=numpy.histogram2d(phi_i,phi_j,bins=12,range=[[-180, 180], [-180, 180]], normed=True)
		(y_n,_)=numpy.histogram(phi_i,bins=x,normed=True)
		(x_n,_)=numpy.histogram(phi_j,bins=y,normed=True)
		(x_n, y_n)=numpy.meshgrid(x_n, y_n)
		mutual=numpy.nansum(H*numpy.log(H/(x_n*y_n)))
		return mutual
	else:
		return 0.0
	
	

def main(dir,total_n_residues,n_iterations):
	print "Welcome Human!"
	c=Client(profile='default')
	dihedral_names = ['chi1','chi2','chi3','phi','psi']
	jobs = []

	for id1 in range(1, total_n_residues):
		for id2 in range(id1, total_n_residues):
			for ic in range(0, n_iterations):
				for i, name1 in enumerate(dihedral_names):
					for name2 in dihedral_names[i:]:	
						# construct the jobs
						if ic == 0:
							job = (name1, id1, name2, id2, False,dir)
						else:
							job = (name1, id1, name2, id2, True,dir)
						jobs.append(job)
	result = dview.map(mutual_information_from_files, *zip(*jobs))
	result.wait()
	all_mutuals = result.get()
	grids = {}			

	for i, job in enumerate(jobs):
		name1 = job[0]
		id1 = job[1]
		name2 = job[2]
		id2 = job[3]
		shuffle = job[4]
		if shuffle:
			grid_name = 'grid_ave_'+ (name1) + '_' +(name2)
			if grid_name not in grids:
				grids[grid_name]=numpy.zeros(total_n_residues*total_n_residues)
			grids[grid_name][((id1-1)*total_n_residues)+id2-1] += all_mutuals[i]
		else:
			grid_name = 'grid_'+ (name1) + '_' +(name2)
			if grid_name not in grids:
				grids[grid_name]=numpy.zeros(total_n_residues*total_n_residues)
			grids[grid_name][((id1-1)*total_n_residues)+id2-1] = all_mutuals[i]

#Notes
#the [((id1-1)*total_n_residues)+id2-1] index compresses a 2d array into a 1d array 
#for example, 1,1 is stored in (0*n+0)
#1,n is stored in (0+n-1)
#n,n is stored in(n*(n-1)+n-1)

	
	final_grid=numpy.zeros(((total_n_residues),(total_n_residues)))
	for i, name1 in enumerate(dihedral_names):
		for name2 in dihedral_names[i:]:	
			grid_name_1 = 'grid_'+ name1+ '_' +name2
			grid_name_2 = 'grid_ave_'+ name1 + '_' +name2
			for id1 in range(1, total_n_residues):
				for id2 in range(id1, total_n_residues):
					if n_iterations > 1:
						final_grid[id1-1][id2-1] += grids[grid_name_1][((id1-1)*total_n_residues)+id2-1] - grids[grid_name_2][((id1-1)*total_n_residues)+id2-1]/n_iterations
						final_grid[id2-1][id1-1] = final_grid[id1-1][id2-1] 
					else:
						final_grid[id1-1][id2-1] += grids[grid_name_1][((id1-1)*total_n_residues)+id2-1]
						final_grid[id2-1][id1-1] = final_grid[id1-1][id2-1] 
	file=open('%s'%dir+'mut-inf-mat.txt', 'w')
	numpy.savetxt(file, final_grid)
	


def parse_commandline():
	parser = optparse.OptionParser()
	parser.add_option('-d', '--directory', dest='dir',default="../", help='directory with chi1,phi,psi files')
	parser.add_option('-t','--total_residues',dest='t',type="int",help='total number of residues')
	parser.add_option('-i','--total_iterations',dest='i',type="int",default=1,help='total iterations')
	(options, args) = parser.parse_args()
	return (options, args)
	
#function main if namespace is main

if __name__ == "__main__":
	(options, args) = parse_commandline()
	main(dir=options.dir, total_n_residues=options.t,n_iterations=options.i)


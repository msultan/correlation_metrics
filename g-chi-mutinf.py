#!/usr/bin/env python
from __future__ import division
import numpy    # for histogram stuff
import optparse
import os
import glob
import random
from collections import defaultdict
from IPython import parallel
import sys
import scipy
import string
from tables import *
import re
import time 

def grid_writer(job,value,final_grid,average_grid):
    name1 = job[0]
    id1 = job[1]
    name2 = job[2]
    id2 = job[3]
    shuffle = job[4]
    #if shuffle is active, then put the result in the average grid 
    if shuffle:
        average_grid[id1-1][id2-1] += value
    	#since only half the data is calculated we need to add it for the other half, the check makes sure we don't doubly add the data
        if id1 != id2:
            average_grid[id2-1][id1-1] += value
    else:
        final_grid[id1-1][id2-1] += value

        if id1 != id2:
            final_grid[id2-1][id1-1] += value

def job_checker(dir,jobs,final_grid,average_grid):
	import csv
	from collections import namedtuple
	jobTuple=namedtuple("jobTuple",["parameter"])
	if os.path.exists(dir+'/temp-list.csv'):
	   print "Previous Checkpoint Found:Pruning job list and writing out the matrix"
	   result_dict = {}
	   for [key, val] in csv.reader(open("temp-list.csv")):
	       result_dict[key] = val
           
               for i,job in enumerate(jobs):
                if result_dict.has_key(str(job)):
                    jobs.pop(i)
		    value = result_dict[str(job)]
		    grid_writer(job,float(value),final_grid,average_grid)	    
		
	return jobs     

def writeToCSV(dir,dict):
	import csv
	w = csv.writer(open(dir+"temp-list.csv", "a"))
	for key in dict:
    	   w.writerow([key,dict[key]])


def create_hd5files_from_xvg(dir,skiprows):
    filename_list=glob.glob('%s/*.xvg'%dir)
    for i,filename in enumerate(filename_list):
    #creating hdf5 files for those that have not been created before. 
        if not (os.path.exists(filename+'.h5')):
            h5file=openFile(filename+'.h5',"w",title="time")
            time_dihedral=h5file.createGroup("/",'time')
            data=numpy.loadtxt(filename,usecols=(1,),skiprows=skiprows)
            grid=h5file.createArray('/time','grid',data)
            h5file.close()

def average(x):
    return float(sum(x)) / len(x)	

def pearson_def(x, y):
    import math
    assert len(x)==len(y) 
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)
    

def mutual_information_from_files(res_name1, res_id1, res_name2, res_id2,           shuffle, dir, skiprows,bin_n, test,numWin):
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
    import numpy    #for histogram stuff
    import os
    import scipy.special as sps
    import scipy.stats as scipy
    from scipy.stats.stats import pearsonr
    import glob
    from scipy.stats.stats import pearsonr
    import tables
    import sys
    result = []
    #creating a list with job in first column and result in 2nd column, this will be useful in restarting jobs
    job=(res_name1, res_id1, res_name2, res_id2,         shuffle, dir, skiprows,bin_n, test,numWin)    
    result.append(job)

    if not test:
            # example file namephiALA1.xvg.all_data_450micro
            #f.root.time_dihedral.grid[1]
        filename_id1 = glob.glob('%s/%s???%d.xvg.h5' %(dir,res_name1,res_id1))
        filename_id2 = glob.glob('%s/%s???%d.xvg.h5' %(dir,res_name2,res_id2))
        if  filename_id1 and filename_id2:
            f = tables.openFile(filename_id1[0])
            dihedral_id1 = numpy.array(f.root.time.grid[:])
                
            f = tables.openFile(filename_id2[0])
            dihedral_id2 = numpy.array(f.root.time.grid[:])
                
            if shuffle:
                dihedral_id1 = numpy.random.permutation(dihedral_id1)
                dihedral_id2 = numpy.random.permutation(dihedral_id2)
        else:
            return result.append(0.0)
            sys.out.flush()
	
    if test:
		mu1 = 0.4
		mu2 = 0.5
		kappa1 = numpy.random.random()
		kappa2 = numpy.random.random()

		dihedral_id1 = 57.2957795*numpy.random.vonmises(mu1,kappa1,100000)
		dihedral_id2 = 57.2957795*numpy.random.vonmises(mu2,kappa2,100000)

#The magic begins here:The mutual information is defined as the 2d  sum of joint probability of x,y multiplied by the log(joint/marginal probabilites probabilities)

    #for windowIndex in range(1,numWin):
    res_id1_range = [-180,180]
    res_id2_range = [-180,180]
    
    (joint_pdf_bin_count,id1_bin_edges,id2_bin_edges) = numpy.histogram2d(dihedral_id1[:len(dihedral_id1)],dihedral_id2,bins=bin_n,range = [res_id1_range, res_id2_range],normed = False)
    dihedral_id1_count = numpy.sum(joint_pdf_bin_count,axis = 0)
    dihedral_id2_count = numpy.sum(joint_pdf_bin_count,axis = 1)
    
    dx = id1_bin_edges[1]-id1_bin_edges[0]
    N = len(dihedral_id1)
    H_x = -numpy.nansum((dihedral_id1_count/N) *numpy.log(dihedral_id1_count/N))+numpy.log(dx)
    H_y = -numpy.nansum((dihedral_id2_count/N) *numpy.log(dihedral_id2_count/N))+numpy.log(dx)
    H_x_y = -numpy.nansum((joint_pdf_bin_count/N) *numpy.log(joint_pdf_bin_count/N))+numpy.log(dx*dx)
    
    mutual=H_x+H_y-H_x_y
#normalizing the results
#based on the fact that I(X,Y) <= min(H(x),H(y)); I am using the geometric mean here.
    #mutual = mutual / np.sqrt(H_x * H_y)
    result.append(mutual)
    
    if not test:
        return result
    	sys.out.flush()

    if test:
    
        p = pearson_def(dihedral_id1,dihedral_id2)
        mutual_a = numpy.log(numpy.sqrt(1/(1-p)))
	H_x_a = numpy.log(2*180*sps.jn(0,kappa1))-kappa1*(sps.jn(1,kappa1)/sps.jn(0,kappa1))
	H_y_a = numpy.log(2*180*sps.jn(0,kappa2))-kappa2*(sps.jn(1,kappa2)/sps.jn(0,kappa2))
	H_x_y_a = (0.5*numpy.log(numpy.linalg.det(numpy.cov(dihedral_id1,dihedral_id2))*(2*numpy.pi*numpy.e)**2))
		##print stuff out
	print "\nTesting on VonMises Distribution with mu1,mu2,kappa1,kappa2",mu1,mu2,kappa1,kappa2,"with",bin_n,"bins\n"
	print "Variable","*"*20, "Calculated","*"*20,"Analytical"
	print "H_x","*"*25,H_x,"*"*20,H_x_a
	print "H_y","*"*25,H_y,"*"*20,H_y_a
	print "H_x_y","*"*24,H_x_y,"*"*20,H_x_y_a
	print "mutual(H_x+H_y-H_x_y)","*"*7,mutual*np.sqrt(H_x*H_y),"*"*15,mutual_a
	print "\nNote that every now and then the result might be a bit off from the analytical solution due to the inherent randomness in the starting dataset.\n"

def main(dir,total_n_residues,n_iterations,skiprows,bin_n, test,numWin,prf):
    olderr = numpy.seterr(all='ignore') 
    #Setting up the parallel stuff 
    client_list = parallel.Client(profile=prf)
    print "Running on:",len(client_list.ids)
    view = client_list.load_balanced_view()
    
    #dihedral names
    dihedral_names = ['chi1','chi2','chi3','chi4','phi','psi']
    
    #final job list. 
    jobs = []
    
    #dict of results
    results_dict = {}
    
    #time_jump for saving results in seconds
    time_jump = 1

    
    if test:
        print "TESTING BEGINNING:"
        mutual_information_from_files('junk',1,'junk',1,True,dir,skiprows,bin_n,True,numWin)
        print "TESTING FINISHED"
        sys.exit()
    else:
        final_grid = numpy.zeros(((total_n_residues),(total_n_residues)))
        average_grid = numpy.zeros(((total_n_residues),(total_n_residues)))
        print "Creating Job List"
        file_name_list=glob.glob('%s/*.h5'%dir)
        print "Found",len(file_name_list),"files"
        #for all possible files found
        for i,file_id1 in enumerate(file_name_list):
        #extract the name of the dihedral 
	    name1=re.findall("chi1|chi2|chi3|chi4|phi|psi",file_name_list[i])[0]
	    #extract the id
            id1=int((re.findall("[A-Z]{3}\d+",file_name_list[i])[0])[3:])
        #limit the file_list to what ever is left so that we dont overcount
            for j,file_id2 in enumerate(file_name_list[i:]):
        #get the name and id for the second residue
		name2=re.findall("chi1|chi2|chi3|chi4|phi|psi",file_name_list[i+j])[0]
                id2=int((re.findall("[A-Z]{3}\d+",file_name_list[j])[0])[3:])
                for ic in range(0,n_iterations):
                    if ic==0:
                        job=(name1, id1, name2, id2, False,dir,skiprows,bin_n,False,numWin)
                    else:
                        job = (name1, id1, name2, id2, True,dir,skiprows,bin_n,False,numWin)
                    jobs.append(job)
    print "Created",len(jobs),"jobs"               
    print "Pruning job list"
    job_checker(dir,jobs,final_grid,average_grid)
    if len(jobs) > 0:
        st=time.time()
        print "Start the jobs with ", len(jobs),"jobs at",st
        result = view.map(mutual_information_from_files, *zip(*jobs),chunksize=10)

        print "Jobs Sent ;Waiting on Results Now"
       
        pending=set(result.msg_ids)
        
        #As long as there are pending jobs,
        #wait a while for the jobs to finish
        while pending:
            try:
                client_list.wait(pending,3600)
            except parallel.TimeoutError:
            	pass
            #   print "Here" 
            finished = pending.difference(client_list.outstanding)
            pending = pending.difference(finished)
            for msg_id in finished:
                ar = client_list.get_result(msg_id)
        	for singleResult in ar.result:
        	    print singleResult
        	    results_dict[(singleResult[0])] = singleResult[1]
            writeToCSV(dir,results_dict)
    
        print "BACKED UP THE DATA at %d"%(int(time.time()))
        
            
        all_mutuals = result.get()
        for i,job in enumerate(jobs):
            grid_writer(job,(all_mutuals[i])[1],final_grid,average_grid)
    
    #the final result is the mutual information minus the average of the noise when the # of iterations is greater than 1
    if n_iterations > 1:
    	final_grid = final_grid - (average_grid/(n_iterations-1)) 
    file = open('%s/'%dir+'mut-inf-mat.txt', 'w')
    numpy.savetxt(file, final_grid)


def parse_commandline():
    import os
    parser = optparse.OptionParser()
    parser.add_option('-d', '--directory', dest='dir',default=os.getcwd(), help='directory with chi1,phi,psi files')
    parser.add_option('-t','--total_residues',dest='t',type="int",help='total       number of residues')
    parser.add_option('-i', '--total_iterations', dest='i', type="int", default=1, help=' total iterations')
    parser.add_option('--test', dest='test', default=0, help='test the code')
    parser.add_option('-s', '--skip_rows', dest='s',type='int', default=12,help='how many rows to skip')
    parser.add_option('-b', '--bins',dest='bin_n', type='int', default=24, help='The number of Bins used to bin data(Default 24). Too few bins can lead to problems')
    parser.add_option('-p','--parallel_profile',dest='prf',default='mpi',help='the parallel profile to use')
    parser.add_option('-w', '--windows',dest='numWin', type='int', default=1, help='Whether or not to use windows in the calculation')
    (options, args) = parser.parse_args()
    return (options, args)

#function main if namespace is main

if __name__ == "__main__":
    (options, args) = parse_commandline()
    #create_hd5files_from_xvg(options.dir,options.s)
    main(dir = options.dir, total_n_residues = options.t,n_iterations = options.i,skiprows = options.s,bin_n=options.bin_n, test = options.test,numWin = options.numWin,prf = options.prf)

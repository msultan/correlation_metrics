#!/usr/bin/env python
from __future__ import division
from scipy import spatial
import numpy as np 
import scipy.special
import numpy as np 
import ann
from scipy.spatial import kdtree

	
def count_dist(N,data,kdist):
	from scipy.spatial.distance import euclidean as euc
	#This routine count how many points in the dataset data have 
	#distance lesser then kdist to the point x.
	_,dim=np.shape(data)
	countMatrix=np.zeros((N,2))
	#for all possible points
	for j in range(0,N):
	#check its neighbors 
		for i in range(0,N):
			#calculate distance in X and Y direction using max norm 
			dTempX = np.linalg.norm(data[i][:dim/2] - data[j][:dim/2],np.inf)
			dTempY = np.linalg.norm(data[i][dim/2:] - data[j][dim/2:],np.inf)
			
			#if distance is less than the distance in hyper space
			if ( dTempX < kdist[j] ):
				countMatrix[j][0] += 1
			
			#same for the other direction
			if ( dTempY < kdist[j] ):
				countMatrix[j][1] += 1	
	#-1 because otherwise we would count the actual pt as well.
	return countMatrix-1

def MI(N,k,kdist,data):

	countMat=count_dist(N,data,kdist)
	'''
	countMat should have N rows and 2 cols, we add the psi values of \
	cols up first and then take the average of the (N,1) vector. 
	'''
	out = - np.mean(np.sum(scipy.special.psi(countMat+1),axis=1))

	#for every pt find its neighbors in x and y regime based on the distance to the
	#kth neighbor
	return scipy.special.psi(k)+scipy.special.psi(N) + out
	

def mutual_nearest_neighbors(N,k,i,j,data):
	#calculate the distance to the kth neighbor for every pt in the dataset 
	#using euclidean distances
	kdist=np.zeros(N)
	kdTree = ann.kd_tree(data)
	idx,dis=kdTree.search(data[:],k+1,eps=0.0)
	#temp
	#print idx,dis
	kdist=dis[:,k]	
	mi=MI(N,k,kdist,data)
	#print mi
	#kdist=np.zeros(N)
	#kdTree=kdtree.KDTree(data)
	#dis,idx=kdTree.query(data[:],k+1,p=np.inf,eps=0.0)
	#print idx,dis
	#kdist=dis[:,k]
	mi=MI(N,k,kdist,data)
	#print mi
	return (i,j,mi)

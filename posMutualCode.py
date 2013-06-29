#!/usr/bin/env python
from __future__ import division
from scipy import spatial
from scipy.spatial import KDTree
import numpy as np 
import scipy.special
import numpy as np 

def kdistANN(N,kN,data):
    ##calculate the distance to the kth neighbor for every pt in the dataset 
    #using euclidean distances
    from scipy.spatial.distance import euclidean as euc
    kdist=np.zeros(N)
    kdTree = spatial.KDTree(data)  
    temp=kdTree.query(data[:],kN+1,p=np.inf)
    temp=temp[0][:]   
    kdist=temp[:,kN]
    return kdist

def count_dist(N,data,kdist,refPt):
    from scipy.spatial.distance import euclidean as euc
    #This routine count how many points in the dataset data have 
    #distance lesser then dist to the point x.
    _,dim=np.shape(data)
    countMatrix=np.zeros((2))  
    for i in range(0,N):
        dTempX=np.linalg.norm(data[i][:dim/2]-refPt[:dim/2],np.inf)
        dTempY=np.linalg.norm(data[i][dim/2:]-refPt[dim/2:],np.inf)
        if(dTempX < kdist):
            countMatrix[0] +=1
        if(dTempY < kdist):
            countMatrix[1] +=1   
    #-1 because otherwise we would count the actual pt as well.
    return countMatrix-1

def MI(N,k,kdist,data):

    out=np.zeros(N)
    for i in range(N):
        countMat=count_dist(N,data,kdist[i],data[i])
        out[i] = -np.sum(scipy.special.psi(countMat+1))

    
    #for every pt find its neighbors in x and y regime based on the distance to the
    #kth neighbor
    return scipy.special.psi(k)+scipy.special.psi(N) + np.mean(out)
    

def mutual_nearest_neighbors(N,k,i,j,data):
    kdist=kdistANN(N,k,data)
    mi=MI(N,k,kdist,data)
    return (i,j,mi)
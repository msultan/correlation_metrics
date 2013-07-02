from __future__ import division
import posMutualCode
from posMutualCode import * 
import numpy as np

def printMessage():
	print"Testing the script for a set of gaussians"
	return

def testPositionMutualInformationCode():
	import matplotlib.pyplot as mlp
	printMessage()
	tN=3
	n_frm=N=600
	k=6
	mean=[[0,1,2,3,4,5,6,7,8,9],[0,1,2,3,4,5,6,7],[0,1,2,3,4,5],[0,1,2,3],[0,1]]
	cov=[[[1,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85],[0.85,1,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85],\
	[0.85,0.85,1.0,0.85,0.85,0.85,0.85,0.85,0.85,0.85],[0.85,0.85,0.85,1,0.85,0.85,0.85,0.85,0.85,0.85],\
	[0.85,0.85,0.85,0.85,1,0.85,0.85,0.85,0.85,0.85],[0.85,0.85,0.85,0.85,0.85,1,0.85,0.85,0.85,0.85],\
	[0.85,0.85,0.85,0.85,0.85,0.85,1.0,0.85,0.85,0.85],[0.85,0.85,0.85,0.85,0.85,0.85,0.85,1.0,0.85,0.85],\
	[0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,1.0,0.85],[0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,1.0]],\
	[[1,0.85,0.85,0.85,0.85,0.85,0.85,0.85],[0.85,1,0.85,0.85,0.85,0.85,0.85,0.85],\
	[0.85,0.85,1.0,0.85,0.85,0.85,0.85,0.85],[0.85,0.85,0.85,1,0.85,0.85,0.85,0.85],\
	[0.85,0.85,0.85,0.85,1,0.85,0.85,0.85],[0.85,0.85,0.85,0.85,0.85,1,0.85,0.85],\
	[0.85,0.85,0.85,0.85,0.85,0.85,1.0,0.85],[0.85,0.85,0.85,0.85,0.85,0.85,0.85,1.0]],
	[[1,0.85,0.85,0.85,0.85,0.85],[0.85,1,0.85,0.85,0.85,0.85],[0.85,0.85,1.0,0.85,0.85,0.85],\
	[0.85,0.85,0.85,1,0.85,0.85],[0.85,0.85,0.85,0.85,1,0.85],[0.85,0.85,0.85,0.85,0.85,1]],
	[[1,0.85,0.85,0.85],[0.85,1,0.85,0.85],[0.85,0.85,1.0,0.85],[0.85,0.85,0.85,1]],
	[[1,0.85],[0.85,1]]]
	dim=[10,8,6,4,2]
	# mean=[1,1,1]
	# cov=[[1,0.85,0.85],[0.85,1,0.85],[0.85,0.85,1]]
	# dim=3
	jobs=zip(mean,cov,dim)
	for mean,cov,dim in jobs:
# 		cov=np.random.random((dim,dim))
#  		cov=np.array(cov)
# 		cov =cov * np.random.
		print cov
		tD=np.zeros(tN)
		print "Testing Position MutualInformation Code for %d iterations for %d\
		dimensional Gaussian" %(tN,dim)
		for t in range(tN):
			#get random data	
			data=np.zeros((N,dim))
			coords=np.random.multivariate_normal(mean,cov,size=n_frm)
			#set it to data
			data[:,:]=coords[:,:]
			#subtract the mean from the data
			#data=data-[np.mean(coords[:,0]),0,0,0,np.mean(coords[:,1]),0]
			#data[:,1]=data[:,3]=0
			#data=np.reshape(data.flatten(),(n_frm,6))
			mutual=mutual_nearest_neighbors(N,k,0,0,data)
			print t,mutual[-1]
			tD[t]=mutual[-1]
		l=[]
		for i in range(int(dim/2)):
			l.append(cov[:int(dim/2)][i][:int(dim/2)])
		detx=np.linalg.det(l)
		detxy=np.linalg.det(cov)
		print "Analytical Solution:",0.5*np.log(detx*detx/detxy)
		print 'Average of %d Trial for %d datapoints is %f with standard deviation \
		%f'%(tN,N, np.mean(tD),np.std(tD))
	#mlp.plot(tD)
	#mlp.axhline(-0.5*np.log(det))
	#mlp.show()

	return

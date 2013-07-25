dihedral-mutinf
===============
This will have everything related to information theoretic measures. Right now I have implemented two ways to 
calculate Mutual Information.
1).Binning data which is implemented in dihedral MutualInformation
2).Calculating Mutual Information using Kraskov's method which is implemented in posMutualCode

Notes:
==============
1). Make sure you test the code before running it. Both scripts contain test cases for calculating mutual information 
for bivarite distributions that have analytical solutions. 
2). The positional code starts to show variation from the analytical solution for d>6 if you don't have enough data or if the value of k/N is not carefully choosen. Look at the Kraskov Paper for more details. 

3). For the posotional code, you will need the scikits.ann python wrapper for the ANN library. The ANN library NEEDS to compiled with the INFINITY NORM. 
4).The scikit ANN wrapper might give some problems due to the numpy.testing module not being found. The following link has more details \n
http://sharpen.engr.colostate.edu/mediawiki/index.php/Installing_ANN

Papers to Cite/Read for Background:
===================================
The code here is based on the following papers. Please cite them if you end up using this code

Christopher L. McClendon, Gregory Friedland, David L. Mobley, Homeira Amirkhani, Matthew P. Jacobson. Journal of Chemical Theory and Computation. 2009, 5 (9), pp 2486-2502.
Christopher L. McClendon, Lan Hua, Gabriela Barriero, Matthew P. Jacobson. Journal of Chemical Theory and Computation. 2012. Article ASAP.DOI: 10.1021/ct300008d


A. Kraskov, H. StÃ¶gbauer, and P. Grassberger, Estimating mutual information Phys. Rev. E 69 (6) 066138 (2004)

Lange OF, GrubmÃ¼ller H. Generalized Correlation for Biomolecular Dynamics. Proteins: Structure, Function, and Bioinformatics 62: 1053-1061 (2006)

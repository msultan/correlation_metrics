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
2). The positional code starts to show variation from the analytical solution for d>6 if you don't have enough data 
points. Using around 10k-20k points should work but it will take some time to calculate the solution. 

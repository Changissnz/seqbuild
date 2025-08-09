# file contains script to initialize a 
# linear congruential generator and 
# measures the generator output using 
# this program's <MultiMetric> class. 
# Measurements are stored in file denoted 
# by variable name F. 

set G = make lcg with 
`(s,m,a,m0)`. 

run G for 100 iter. 
set V = run G for `200` iter. 

set M = make multimetric with V. 

set X = run M with ngram=5.
set F = open file with `foldername/filename`. 
write X to F. 

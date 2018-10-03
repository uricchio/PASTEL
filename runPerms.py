from PolySel import statTests
import sys

t0 = float(sys.argv[2])
t1 = float(sys.argv[3])

height = statTests.Pheno(numbins=100,infoFilt=0.8)

height.getSize(sys.argv[1])

height.getData(sys.argv[1],t0,t1)

# this line prints the value of S_beta 
# to print the beta-DAF curve, change to 'pr=True'
height.runIter(perm=False,pr=False)

# This loop runs the permutations
for i in xrange(2000):
    height.runIter(perm=True,pr=False)


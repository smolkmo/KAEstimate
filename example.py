from kaestimate import *

start=time.time()

####################################
#BEGIN PARAMETERS
####################################

#Adapt these to your needs

#Database length
m=1000

#Query length
n=100

#Alphabet
alph="ACGT"

#Scoring scheme
match=10
mismatch=-15
gapopen=-20
gapextend=-20


#Computational effort
threads=4
maxtime=60*10 #60*60*13

####################################
#END PARAMETERS
####################################

#Gapped
distr=getBestScoreDistribution(m,n,alph,lambda a,b: getLocalAlignmentScore(a,b,match,mismatch,gapopen,gapextend),threads,maxtime)

#Use this instead for ungapped
#distr=getBestScoreDistribution(m,n,alph,lambda a,b: getLocalUngappedAlignmentScore(a,b,match,mismatch),threads,maxtime)

report=open("report.txt","w")
report.write("Karlin-Altschul Parameter Estimation Report\n")
report.write("===========================================\n")
report.write("Run stats:\n")
report.write("\tm=%d (database); n=%d (query); alph=%s\n"%(m,n,alph))
report.write("\tmatch=%d; mismatch=%d; gapopen=%d; gapextend=%d\n"%(match,mismatch,gapopen,gapextend))
report.write("\titers=%d; threads=%d; maxtime=%d; wallclock=%d\n"%(iters(distr),threads,maxtime,time.time()-start))
report.write("===========================================\n")
report.write("Estimated score p-Values and associated KA parameter estimation:\n")

last_i=None
for i in sorted(distr):
    if last_i==None:
        report.write("\tp(%d)=%f; lambda=????????????; Kappa=????????????\n"%(i,p(distr,i)))
    else:
        try:
            l=getLambdaInternal(m*n,last_i,p(distr,last_i),i,p(distr,i))
            K=getKappaInternal(m*n,last_i,p(distr,last_i),i,p(distr,i))
            report.write("\tp(%d)=%f; lambda=%.10f; Kappa=%.10f (Fit avg. dist: %.3f)\n"%(i,p(distr,i),l,K,getFitDist(distr,l,K,m,n)))
        except ValueError:
            report.write("\tp(%d)=%f; lambda=????????????; Kappa=????????????\n"%(i,p(distr,i)))
    last_i=i

from scipy.optimize import minimize
import numpy as np
def writeBestFit(maxp=1.0):
    report.write("===========================================\n")
    report.write("Best fitting KA parameters (p<=%f):\n"%maxp)
    all_lk=[]
    for i in sorted(distr):
        for j in sorted(distr):
            try:
                l=getLambdaInternal(m*n,last_i,p(distr,last_i),i,p(distr,i))
                K=getKappaInternal(m*n,last_i,p(distr,last_i),i,p(distr,i))
                all_lk.append((l,K,i,j))
            except:
                pass
       
    best_dist=float("inf")

    for l, K, a, b in all_lk:
        try:
            dist=getFitDist(distr,l,K,m,n,maxp)
            if dist<best_dist:
                best_dist=dist
                best_l,best_K=(l,K)
                best_a,best_b=(a,b)
        except:
            return

    report.write("\tDirect:      lambda=%.10f, Kappa=%.10f (Fit avg. dist: %.10f)\n"%(best_l,best_K,getFitDist(distr,best_l,best_K,m,n,maxp)))

    def distance(x):
        return getFitDist(distr,x[0],x[1],m,n,maxp)

    x0=np.array([0.0,0.0])

    res=minimize(distance,x0,method="nelder-mead",options={"disp":True,"xtol":1e-8})

    report.write("\tNelder-Mead: lambda=%.10f, Kappa=%.10f (Fit avg. dist: %.10f)\n"%(res.x[0],res.x[1],getFitDist(distr,res.x[0],res.x[1],m,n,maxp)))


writeBestFit(1.0)
writeBestFit(0.75)
writeBestFit(0.5)
writeBestFit(0.25)
writeBestFit(0.05)

report.flush()
report.close()

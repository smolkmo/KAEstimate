from Bio import pairwise2
import random
import argparse
from multiprocessing import Process, Queue
import time
from math import log,exp

def getSingleSearchBestScore(m,n,alphabet,getAlignmentScoreFn):
    db="".join([random.choice(alphabet) for i in range(m)])
    q="".join([random.choice(alphabet) for i in range(n)])
    max_score=-float("inf")
    for i in range(m):
        db_sub=db[i:i+n]
        if len(db_sub)<n:
            continue

        score=getAlignmentScoreFn(q,db_sub)
        max_score=max(max_score,score)
    return max_score

def worker(q,m,n,alphabet,getAlignmentScoreFn):
    while True:
        score=getSingleSearchBestScore(m,n,alphabet,getAlignmentScoreFn)
        q.put(score)


def getBestScoreDistribution(m,n,alphabet,getAlignmentScoreFn,threads=1,maxtime=60):
    q=Queue()
    workers=[]
    for i in range(threads):
        p=Process(target=worker,args=(q,m,n,alphabet,getAlignmentScoreFn))
        p.start()
        workers.append(p)

    start=time.time()
    best_scores={}
    iters=0
    while time.time()-start < maxtime:
        score=q.get()
        if not score in best_scores:
            best_scores[score]=0
        best_scores[score]+=1
        iters+=1
        print("%d iters, %d%% of time expired"%(iters, (time.time()-start)/(maxtime/100.0) ))
        print(best_scores)
        print("")

    for p in workers:
        p.terminate()

    return best_scores

def getLocalAlignmentScore(a,b,match,mismatch,gap_open,gap_extend):
    return pairwise2.align.localms(a,b,match,mismatch,gap_open,gap_extend,score_only=True)

def p(distr,observed_score):
    sum_iters=0
    for score in distr:
        sum_iters+=distr[score]

    geq_iters=0
    for score in distr:
        if score >= observed_score:
            geq_iters+=distr[score]

    return float(geq_iters)/float(sum_iters)

def minscore(distr):
    return int(min(distr.keys()))

def maxscore(distr):
    return int(max(distr.keys()))

def iters(distr):
    sum_iters=0
    for score in distr:
        sum_iters+=distr[score]
    return sum_iters    

def getLambdaInternal(n,Sa,pa,Sb,pb):
    return -log(log(-pa + 1)/log(-pb + 1))/(Sa - Sb)

def getKappaInternal(n,Sa,pa,Sb,pb):
    return -exp(-Sa*log(log(-pa + 1)/log(-pb + 1))/(Sa - Sb))*log(-pa + 1)/n

def getKAEValue(S,m,n,K,lambd):
    return K*m*n*exp((-lambd)*S)

def getKAPValue(S,m,n,K,lambd):
    return 1-exp(-getKAEValue(S,m,n,K,lambd))

def getFitDist(distr,l,K,m,n,maxp=1.0):
    dist=0
    considered=0

    for i in distr:
        p_act=p(distr,i)
        p_ka=getKAPValue(i,m,n,K,l)

        if p_act <= maxp:
            considered+=1
            dist+=(p_act-p_ka)**2
    return dist/float(considered)
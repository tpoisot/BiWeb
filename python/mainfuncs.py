##########
# bipy
# main file
##########
# A set of python functions to work on bipartite interaction networks
# Timothee Poisot - Universite Montpellier 2
# First version : august 12, 2011
# Contact : timothee.poisot@uqar.ca
##########
import scipy as sp
import numpy as np
import random
##########

def meanperf(W,novoid=False):
    # Mean performance of all TL species
    perf = []
    fit = 0
    for fit in W:
        tspe = []
        for tfit in fit:
            if novoid:
                if tfit>0:
                    tspe.append(tfit)
            else:
                tspe.append(tfit)
        perf.append(np.mean(tspe))
    return perf


def surp(V):
    # Given a vector, returns the
    # probability of each element being a
    # surprise in the sens of Shannon's
    # self-information
    #
    # Elements are rounded to 1 decimal place
    # to get meaningful results
    Nbits = len(V)
    sur = []
    for v in V:
        tsum = 0
        for w in V:
            if round(w,2) == round(v,2):
                tsum += 1
        tsum = float(tsum)/Nbits
        sur.append(tsum)
    return sur


def eH(V):
    # Exponent of Shannon's entropy
    tV = []
    for i in range(len(V)):
        tV.append(round(V[i],2))
    p = surp(tV)
    # Only unique elements
    uV = [tV[0]]
    uP = [p[0]]
    for i in range(1,len(tV)):
        isUnique = 0
        for u in range(0,len(uV)):
            if tV[i] == uV[u]:
                isUnique += 1
                break
        if isUnique == 0:
            uV.append(tV[i])
            uP.append(p[i])
    # If all values are equal, no information
    if len(uV) == 1:
        return float(0)
    # Loop
    pLNp = 0
    for i in range(0,len(uV)):
        if uP[i] == 0:
            pLNp += 0
        else:
            pLNp += uP[i]*np.log(uP[i])
    Hprime = 0 - pLNp - np.log(len(uV))
    return np.exp(Hprime)


def qrange(V):
    # Difference between maximal and minimal values
    # of a vector
    return max(V)-min(V)


def d2h(V,bin='sturgis',prop=True):
    # Returns the values needed to draw an histogram in PyX
    #
    # bin can take the values 'sturgis' or 'rice'
    nInd = len(V)
    ## Number of bins in the histogram
    if bin == 'sturgis':
        nbin = np.ceil(1+np.log2(len(V)))
    if bin =='rice':
        nbin = np.ceil(2*pow(len(V),1.0/3))
    ## General parameters
    bottom = round(min(V),1)
    top = round(max(V),1)
    ran = top - bottom
    span = ran / float(nbin)
    ## Do the data
    xs = []
    ys = []
    # Built the bins
    for i in range(int(nbin)):
        brange = bottom + (i*span)
        trange = bottom + (i+1)*span
        rmean = (brange+trange)/2
        tcount = 0
        for j in range(len(V)):
            if (brange <= V[j]) & (V[j] < trange):
                tcount += 1
        if prop:
            ys.append(tcount/float(nInd))
        else:
            ys.append(tcount)
        xs.append(rmean)
    # And finally...
    return zip(xs,ys)


def count(V):
    # Counts the number of elements with a given integer value
    # (mostly) useful to generate degree plots
    cnt = []
    m = min(V)
    M = max(V)
    r = range(0,(M+1))
    V = sorted(V)
    for ra in r:
        tcnt = 0
        for i in range(0,len(V)):
            if V[i] == ra:
                tcnt += 1
        cnt.append(tcnt)
    return zip(r,cnt)


def spread(V,m,M):
    mi = min(V)
    for i in range(len(V)):
        V[i] = float(V[i])-mi
    ma = max(V)
    for i in range(len(V)):
        V[i] = float(V[i])/ma
        V[i] = V[i]* (M - m)
        V[i] = V[i] + m
    return V



def uniquify(seq):
    """
    Return an array of the unique elements of an array
    """
    set = {}
    map(set.__setitem__, seq, []) 
    return sorted(set.keys())

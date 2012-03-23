import scipy as sp
import numpy as np

from ..mainfuncs import *
from ..base import *
from ..spe import *

def null_1(tem):
    tint = np.prod(tem.shape)
    zeroes = tint-np.sum(tem)
    hasNull = True
    while hasNull:
        mat = np.concatenate((np.ones(np.sum(tem)),np.zeros(zeroes)),1)
        np.random.shuffle(mat)
        mat = mat.reshape((len(tem),len(tem[0])))
        hasNull = (np.sum(mat.sum(axis=0)==0) + np.sum(mat.sum(axis=1)==0)) > 0
    return mat


def null_2(tem):
    Ri = np.sum(tem,axis=1) / float(len(tem[0]))
    Ci = np.sum(tem,axis=0) / float(len(tem))
    hasNull = True
    while hasNull:
        mat = np.zeros(np.prod(tem.shape)).reshape((len(tem),len(tem[0])))
        for row in xrange(len(tem)):
                Pi = (Ri[row]+Ci)/float(2)
                Ro = np.random.binomial(1,Pi)
                mat[row] = Ro
        hasNull = (np.sum(mat.sum(axis=0)==0) + np.sum(mat.sum(axis=1)==0)) > 0
    return mat


def null_3row(tem):
    Ri = np.sum(tem,axis=1) / float(len(tem[0]))
    hasNull = True
    while hasNull:
        mat = np.zeros(np.prod(tem.shape)).reshape((len(tem),len(tem[0])))
        for row in range(len(tem)):
                Ro = np.random.binomial(1,np.ones(len(tem[0]))*Ri[row])
                mat[row] = Ro
        hasNull = (np.sum(mat.sum(axis=0)==0) + np.sum(mat.sum(axis=1)==0)) > 0
    return mat


def null_3col(tem):
    Ci = np.sum(tem,axis=0) / float(len(tem))
    hasNull = True
    while hasNull:
        mat = np.zeros(np.prod(tem.shape)).reshape((len(tem),len(tem[0]))).T
        for col in range(len(tem[0])):
                Co = np.random.binomial(1,np.ones(len(tem))*Ci[col])
                mat[col] = Co
        mat = mat.T
        hasNull = (np.sum(mat.sum(axis=0)==0) + np.sum(mat.sum(axis=1)==0)) > 0
    return mat


def nullModel(template,fun=null_1,replicates=1000,minconn=1.5):
    ad = np.copy(template.adjacency)
    nullWebs = []
    # Check the level of constraints
    numInt = np.sum(ad)
    minInt = min(ad.shape)
    rOrig = False
    if numInt < minInt*minconn:
        print 'WARNING: The connectance is low, you network is highly constrained'
        print '         The original network will be returned'
        rOrig = True
    # Analyses
    for i in range(replicates):
        if rOrig:
            nullWebs.append(ad)
        else:
            nullWebs.append(fun(ad))
    return nullWebs
import scipy as sp
import numpy as np
from ..mainfuncs import *
from ..spe import *
from ..base import *
from scipy.weave import inline

def nestadj(aW):
    return adjacency(sortbydegree(aW))

def compareones(w1,w2,tn):
    #return (100 * np.sum( (w1 == np.ones(len(w1))) * (w2 == np.ones(len(w2)))) )/tn
    return (100 * np.sum( (w1+w2) == 2 ))/tn

def compareones_c(w1,w2,tn):
    """
    A C implementation of the compareones
    """
    nw1 = np.int_(np.copy(w1))
    nw2 = np.int_(np.copy(w2))
    code = """
    int s;
    s = 0;
    for(int i = 0; i < n; i++)
    {
        if((nw1[i] == 1)&(nw2[i] == 1))
         {
            s += 100;
         }
    }
    return_val = s;
    """
    n = len(w2)
    res = inline(code, ['nw1','nw2','n'], headers = ['<math.h>'], compiler = 'gcc')
    return res / float(tn)

def getNpaired(W,strict,use_c):
    # Required for NODF calculation
    # Get the N paired value of a web
    Npaired = []
    gen = generality(W)
    if use_c:
        COfunc = compareones_c
    else:
        COfunc = compareones
    for i in xrange(0,(len(W)-1)):
        for j in xrange((i+1),len(W)):
            if strict:
                if gen[i] > gen[j]:
                    Npaired.append(COfunc(W[i],W[j],gen[j]))
                else :
                    Npaired.append(0)
            else:
                if gen[i] >= gen[j]:
                    Npaired.append(COfunc(W[i],W[j],gen[j]))
                else :
                    Npaired.append(0)
    return Npaired


def nodf(aW,strict=True,use_c=True):
    """
    The strict boolean tells if the condition for overlap needs to be
    strictly enforced or not. If srict = False, this can somehow increase the NODF
    values.
    """
    if (len(aW[0])==1)|(len(aW)==1):
        return [0,0,0]
    W = nestadj(aW)
    NProw = getNpaired(W,strict,use_c)
    NPcol = getNpaired(W.T,strict,use_c)
    # Output the NODF value
    ColCor = (len(W[0])*(len(W[0])-1))/2
    RowCor = (len(W)*(len(W)-1))/2
    NrowSum = 0
    for p in xrange(0,len(NProw)):
        NrowSum += NProw[p]
    NcolSum = 0
    for p in xrange(0,len(NPcol)):
        NcolSum += NPcol[p]
    WholeNest = np.round((NcolSum+NrowSum)/(ColCor + RowCor),2)
    ColNest = np.round(NcolSum/ColCor,2)
    RowNest = np.round(NrowSum/RowCor,2)
    NEST = [WholeNest, ColNest, RowNest]
    return NEST

import scipy as sp
import numpy as np

def adjacency(W):
    """
    Return the adjacency matrix
    """
    adj = np.copy(W)
    adj[adj>0] = 1
    return adj


def websize(W):
    """
    Return the size of a matrix
    """
    return np.prod(np.shape(W))


def connectance(W):
    """
    Return the connectance as L/S^2
    """
    ad = adjacency(W)
    ws = websize(W)
    return float(np.sum(ad))/ws


def generality(W):
    ad = adjacency(W)
    gen = np.sum(ad,axis=1)
    return gen


def vulnerability(W):
    ad = adjacency(W)
    gen = np.sum(ad,axis=0)
    return gen


def rank(V):
    # Returns the rank of a vector
    # with no ties
    rn = np.zeros(len(V),dtype=np.int32)
    crnk = 0
    while crnk < len(V):
        for j in xrange(0,len(V)):
            cMax = np.max(V)
            if V[j] == cMax:
                rn[j] = crnk
                crnk += 1
                V[j] = np.min(V)-1
                break
    return rn


def sortbydegree(W):
    # Sort a matrix by degree
    # Better for visualization
    # Required for nestednes
    if hasattr(W,'connectance'):
        g = W.generality
        v = W.vulnerability
        upsp = W.upsp
        losp = W.losp
        web = W.web
        ## VOID VECTORS FOR THE SORTED SPECIES NAMES
        oTnames = W.upnames
        oBnames = W.lonames
        vTnames = np.copy(oTnames)
        vBnames = np.copy(oBnames)
    else:
        g = generality(W)
        v = vulnerability(W)
        upsp = len(W)
        losp = len(W[0])
        web = W
    ## Step 1 : sort TLO
    rG = rank(g)
    nW = np.zeros((upsp,losp))
    for ro in range(0,upsp):
        nW[rG[ro]] = web[ro]
        if hasattr(W,'connectance'):
            vTnames[rG[ro]] = oTnames[ro]
    ## Step 2 : sort BLO
    nW = nW.T
    dW = np.zeros((upsp,losp)).T
    rG = rank(v)
    for ro in range(0,losp):
        dW[rG[ro]] = nW[ro]
        if hasattr(W,'connectance'):
            vBnames[rG[ro]] = oBnames[ro]
    Fweb = np.copy(dW.T)
    # This is an horrible horrible solution
    if hasattr(W,'connectance'):
        Fweb = [Fweb,vTnames,vBnames]
    return Fweb


def fixmat(aW):
    import numpy as np
    W = np.copy(aW)
    # Fix a matrix so that there are no empty row
    # or empty columns, issues a message is some rows
    # were removed
    OrigSize = websize(W)
    g = generality(W)
    v = vulnerability(W)
    emptyRows = 0
    emptyCols = 0
    for i in range(len(g)):
        if g[i]==0:
            emptyRows += 1
    for j in range(len(v)):
        if v[j]==0:
            emptyCols += 1
    if (emptyRows+emptyCols)==0:
        return W
    else :
        # Create a new matrix with the correct dimensions
        nW = np.zeros(((len(W)-emptyRows),(len(W[0])-emptyCols)),float)
        # For each row and each column
        # copy the correct values in the new matrix
        cRow = 0
        for i in range(len(g)):
            # If the species i is not interacting, we can skip
            if g[i] == 0:
                continue
            else:
                cCol = 0
                # Else we go throug the species j
                for j in range(len(v)):
                    if v[j] > 0:
                        nW[cRow][cCol] = W[i][j]
                        cCol += 1
                cRow += 1
        # Finally...
        return nW


def readweb(fname):
    # Read a web from a text matrix with top trophic level organisms as rows
    data = np.loadtxt(fname)
    return fixmat(data)


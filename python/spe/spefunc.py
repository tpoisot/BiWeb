import scipy as sp
import numpy as np

from ..mainfuncs import *
from ..base import *

def sp_pdi(f):
    fit = np.copy(f)
    tspe = 0
    fit = sorted(fit,reverse=True)
    ma = max(fit)
    for i in xrange(1,len(fit)):
        tspe += (((fit[0]/float(ma))-(fit[i]/float(ma)))/(len(fit)-1))
    return tspe


def specificity(W):
    """
    Returns the Paired Differences Index measure of specificity
    Normalized between 0 and 1
    """
    spe = []
    for fit in W:
        spe.append(round(sp_pdi(fit),3))
    return spe


def ssi(W):
    # Measures the Species Specialization Index
    # usefull to discriminate among highly specialized organisms
    #	Julliard, R, et al. (2006) Ecol Lett 9(11) 1237-1244 10.1111/j.1461-0248.2006.00977.x
    spe = []
    n = len(W[0]) # Number of resources
    normfac = n * np.sqrt((n-1)/float(n))
    for fit in W:
        # For each species in the web
        tspe = 0
        Pbar = np.mean(fit)
        tnfac = normfac * Pbar
        for sp in fit:
            tspe += (sp-Pbar) ** 2
        spe.append(round(np.sqrt(tspe)/tnfac,3))
    return spe


def IR(W):
    # Correlation matrix of responses
    rho = np.corrcoef(W.web)
    # Standard deviations of responses
    sig = []
    for fit in W.web:
        sig.append(round(np.std(fit),4))
    # Correction factor
    cfac = W.upsp*(W.upsp-1)
    # Measures
    tR = 0
    tI = 0
    for i in range(W.upsp-1):
        for j in range((i+1),W.upsp):
            tR += (sig[i]-sig[j]) ** 2
            tI += sig[i]*sig[j]*(1-rho[i][j])
    R = tR / (2 * cfac)
    I = tI / cfac

    VarW = np.std(W.web)*float(np.std(W.web))

    if VarW == 0:
        R = 0
        I = 0
    else:
        R = R / float(VarW)
        I = I / float(VarW)
    return [R,I]


def RR(W):
    ge = generality(W)
    R = len(W[0])
    rr = []
    for i in range(len(ge)):
        rr.append((R-ge[i])/float(R-1))
    return rr


def HighLink(W):
    hl = []
    for i in range(len(W)):
        hl.append(max(W[i]))
    return hl

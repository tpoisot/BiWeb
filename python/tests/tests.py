from ..null import *
from ..nes import *
from ..mod import *
import scipy.stats as spp


def PValToText(pv):
    ptxt = ' --- '
    if pv < 0.05:
        ptxt = '  *  '
    if pv < 0.01:
        ptxt = '  ** '
    if pv < 0.001:
        ptxt= ' *** '
    if pv < 0.0001:
        ptxt = '**** '
    if pv < 0.00001:
        ptxt = '*****'
    return ptxt

class test:
    def __init__(self,web):
        """
        Initialize a test class and perform the tests
        """
        self.web = web
        self.nulls = []
        self.devnest = []
        self.devnest_lo = []
        self.devnest_up = []
        self.devqb = []
        self.devqr = []
        self.model = null_1
        self.replicates = 100
    def donulls(self,model=null_1,replicates=100):
        self.model = model
        self.replicates = replicates
        self.nulls = nullModel(self.web,self.model,self.replicates)
    def nestedness(self):
        if len(self.nulls) == 0:
            self.donulls()
        d_nest = getDevNest(self.web,self.nulls,self.web.use_c)
        self.devnest = d_nest[0]
        self.devnest_lo = d_nest[1]
        self.devnest_up = d_nest[2]
    def modularity(self,repl):
        if len(self.nulls) == 0:
            self.donulls()
            ## test if the bipartite object has modules
        if not self.web.modules.done:
            self.web.modules.detect(reps=repl,use_c=self.web.use_c)
            ##
        d_mod = getDevMod(self.web,self.nulls,repl,self.web.use_c)
        self.devqr = d_mod[0]
        self.devqb = d_mod[1]
    def __str__(self):
        out = "Stat\tN0\tN'\tp\tIC-\tIC+\n"
        out +=  "---------------------------------------------\n"
        if len(self.devnest) > 0:
            out += " NODF\t"+str(round(self.devnest[0],2)).zfill(4)+"\t"+str(round(self.devnest[2],2)).zfill(4)+"\t"+PValToText(self.devnest[1])+"\t"+str(round(self.devnest[3],2)).zfill(4)+"\t"+str(round(self.devnest[4],2)).zfill(4)+"\n"
        if len(self.devnest_lo) > 0:
            out += "bNODF\t"+str(round(self.devnest_lo[0],2)).zfill(4)+"\t"+str(round(self.devnest_lo[2],2)).zfill(4)+"\t"+PValToText(self.devnest_lo[1])+"\t"+str(round(self.devnest_lo[3],2)).zfill(4)+"\t"+str(round(self.devnest_lo[4],2)).zfill(4)+"\n"
        if len(self.devnest_up) > 0:
            out += "tNODF\t"+str(round(self.devnest_up[0],2)).zfill(4)+"\t"+str(round(self.devnest_up[2],2)).zfill(4)+"\t"+PValToText(self.devnest_up[1])+"\t"+str(round(self.devnest_up[3],2)).zfill(4)+"\t"+str(round(self.devnest_up[4],2)).zfill(4)+"\n"
        if len(self.devqr) > 0:
            out += " QR    \t"+str(round(self.devqr[0],2)).zfill(4)+"\t"+str(round(self.devqr[2],2)).zfill(4)+"\t"+PValToText(self.devqr[1])+"\t"+str(round(self.devqr[3],2)).zfill(4)+"\t"+str(round(self.devqr[4],2)).zfill(4)+"\n"
        if len(self.devqb) > 0:
            out += " QB    \t"+str(round(self.devqb[0],2)).zfill(4)+"\t"+str(round(self.devqb[2],2)).zfill(4)+"\t"+PValToText(self.devqb[1])+"\t"+str(round(self.devqb[3],2)).zfill(4)+"\t"+str(round(self.devqb[4],2)).zfill(4)+"\n"
        return out

def gMIC(distrib):
    """
    Bayesian estimates of mean, and standard deviation
    """
    estimates = spp.bayes_mvs(distrib,alpha=0.95)[0]
    return[estimates[0],estimates[1][0],estimates[1][1]]

def getDevNest(w,list,use_c):
    expect = []
    expect_up = []
    expect_lo = []
    for i in list:
        Nodf = nodf(i,strict=w.nodf_strict,use_c=use_c)
        expect.append(Nodf[0])
        expect_up.append(Nodf[2])
        expect_lo.append(Nodf[1])
    testRes = spp.ttest_1samp(expect, w.nodf)
    testRes_up = spp.ttest_1samp(expect_up, w.nodf_up)
    testRes_lo = spp.ttest_1samp(expect_lo, w.nodf_low)
    OUT = [w.nodf,testRes[1]]
    OUT_up = [w.nodf_up,testRes_up[1]]
    OUT_lo = [w.nodf_low,testRes_lo[1]]
    est = gMIC(expect)
    est_lo = gMIC(expect_lo)
    est_up = gMIC(expect_up)
    for est_par in est:
        OUT.append(est_par)
    for est_par in est_lo:
        OUT_lo.append(est_par)
    for est_par in est_up:
        OUT_up.append(est_par)
    return [OUT, OUT_lo, OUT_up]

def getDevMod(w,nulls,rep,use_c):
    """
    Get the deviation from random expectation of modularity. Optimized so that
    the null webs are gone through only one time. Retunrs two arrays, one for
    Qr, the other for Qb.
    """
    m = [w.modules.Q,w.modules.N,w.modules.up_modules,w.modules.low_modules]
    Qbsim = []
    Qrsim = []
    wQr = Qr(w.web,m)
    wQb = w.modules.Q
    for c_null in nulls:
        c_mod = findModules(c_null, use_c = use_c)
        Qrsim.append(Qr(c_null, c_mod))
        Qbsim.append(c_mod[0])
    testResB = spp.ttest_1samp(Qbsim, wQb)
    testResR = spp.ttest_1samp(Qrsim, wQr)
    OUT_r = [wQr,testResR[1]]
    OUT_b = [wQb,testResB[1]]
    est_r = gMIC(Qrsim)
    est_b = gMIC(Qbsim)
    for est_par in est_r:
        OUT_r.append(est_par)
    for est_par in est_b:
        OUT_b.append(est_par)
    return [OUT_r,OUT_b]
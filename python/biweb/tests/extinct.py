import numpy as np
import copy
from ..mainfuncs import *
from numpy.random import shuffle

class ext_seq:
    def __init__(self,name='random'):
        self.name = name
        self.result = []
        self.score = 0
        self.removed = []
        self.survived = []

class robustness:
    def __init__(self,web):
        self.web = web
        self.random = ext_seq('random')
        self.stog = ext_seq('stog')
        self.gtos = ext_seq('gtos')
    def do_random(self,lower=True,reps=100):
        temp = extinctRobustness(self.web,self.random.name,lower,reps)
        self.random.result = temp
        sc_temp = extinctionScore(self.random.result)
        self.random.score = sc_temp[0]
        self.random.removed = sc_temp[1]
        self.random.survived = sc_temp[2]
    def do_stog(self,lower=True,reps=5):
        temp = extinctRobustness(self.web,self.stog.name,lower,reps)
        self.stog.result = temp
        sc_temp = extinctionScore(self.stog.result)
        self.stog.score = sc_temp[0]
        self.stog.removed = sc_temp[1]
        self.stog.survived = sc_temp[2]
    def do_gtos(self,lower=True,reps=5):
        temp = extinctRobustness(self.web,self.gtos.name,lower,reps)
        self.gtos.result = temp
        sc_temp = extinctionScore(self.gtos.result)
        self.gtos.score = sc_temp[0]
        self.gtos.removed = sc_temp[1]
        self.gtos.survived = sc_temp[2]
    def __str__(self):
        s = 'RND: '+str(self.random.score)+'\n'
        s+= 'S>G: '+str(self.stog.score)+'\n'
        s+= 'G>S: '+str(self.gtos.score)
        return s
    def plot(self):
        try:
            import pyx
        except ImportError:
            print 'PyX is not installed on your system'
            print 'You can get it at http://pyx.sourceforge.net/'
            print 'Install PyX before using this function -- STOP'
            return 0
        g = pyx.graph.graphxy(width=8,
                          x = pyx.graph.axis.linear(title='Proportion removed'),
                          y = pyx.graph.axis.linear(title='Generality'),
                          key = pyx.graph.key.key(pos="bl", dist=0.1))
        g.plot([
            pyx.graph.data.values(x=self.stog.removed, y=self.stog.survived, title='Increasing'),
            pyx.graph.data.values(x=self.random.removed, y=self.random.survived, title='Random'),
            pyx.graph.data.values(x=self.gtos.removed, y=self.gtos.survived, title='Decreasing'),
            ],[pyx.graph.style.line()])
        g.writePDFfile(self.web.name+"_ext_degree")
        return 1


def extinctRobustness(w,method='random',removelower=True,nreps=100):
    ######
    #
    # method
    # Must be one of
    #	random : random removal
    #	stog : specialists to generalists
    #	gtos : generalists to specialists
    #
    # removelower
    #	tells if the upper or lower level is removed
    #
    # nreps
    #	Number of extinction sequences to perform
    #
    ######

    # Keep the data somewhere
    # and make sure that the extreme points are accounted for
    NSP_REMV = [0,1]
    NSP_SURV = [1,0]

    # Actual simulations
    currentRep = 0

    # Loop for the extinction sequence
    while currentRep < nreps:
        # Number of links for stog/gtos sequences
        if removelower:
            Links = w.vulnerability
            Range = np.array(range(0,w.losp))
            fTop = False
        else:
            Links = w.generality
            Range = np.array(range(0,w.upsp))
            fTop = True

        inds = np.argsort(np.copy(Links))
        Range = np.take(Range, inds)

        if method == 'gtos':
                Range = Range[::-1]

        # Increase the sequence counter
        currentRep += 1
        # Randomize if needed
        if method == 'random':
            shuffle(Range)
        # Begin extinctions
        mat = np.copy(w.adjacency)
        for i in range(0,(len(Range)-1)):
            if fTop:
                for j in range(0,w.losp):
                    mat[Range[i],j] = 0
            else:
                for j in range(0,w.upsp):
                    mat[j,Range[i]] = 0
            if removelower:
                NSP_REMV.append((i+1)/float(w.losp))
                NSP_SURV.append(np.sum(mat.sum(axis=1)>0)/float(w.upsp))
            else:
                NSP_REMV.append((i+1)/float(w.upsp))
                NSP_SURV.append(np.sum(mat.sum(axis=0)>0)/float(w.upsp))

    # Finally...
    return zip(NSP_REMV,NSP_SURV)


def int_rect(x,y):
    # Numerical integration using rectangles method
    AUC = 0
    LastLowerBound = 0
    LastUpperBound = 0
    for i in range(len(x)):
        # Find span and bounds
        if LastUpperBound == LastLowerBound:
            Span = x[(i+1)]-x[i]
            LastUpperBound = x[i]+Span/float(2)
            LastLowerBound = x[i]-Span/float(2)
        else:
            LastLowerBound = LastUpperBound
            HSpan = (x[i]-LastLowerBound)
            Span = HSpan  *2
            LastUpperBound = LastLowerBound + Span
        AUC += y[i]*Span
    return AUC

def extinctionScore(ext,integrator=int_rect):
    # Gives the area under the curve for an extinction robustness analysis
    # Step 1 : aggregate
    NRem = list(zip(*ext)[0])
    NSur = list(zip(*ext)[1])
    #
    MSur = []
    LRem = uniquify(NRem)
    #
    for lrem in LRem:
        tsum = 0
        tcnt = 0
        for i in range(len(NSur)):
            if NRem[i] == lrem:
                tsum += NSur[i]
                tcnt += 1
        MSur.append(tsum/float(tcnt))
    # Step 2 : perform integration
    AUC = integrator(LRem,MSur)
    return [AUC,LRem,MSur]


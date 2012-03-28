from .nes import *
from .mod import *
from .graphs import *
from .null import *
from .tests import *

import pickle
import tempfile
import os
import numpy as np
import networkx as nx

import os.path

class bipartite:
    ## This class defines a bipartite object with all structural infos
    def __init__ (self,web,t=False,nodf_strict=True,use_c=True):
        self.use_c = use_c
        # Read the matrix
        if t:
            web = web.T
        self.web = fixmat(web)
        # General infos
        self.upsp = len(self.web)
        self.losp = len(self.web[0])
        self.size = self.upsp * self.losp
        # Connectance
        self.adjacency = adjacency(self.web)
        self.nlink = np.sum(self.adjacency)
        self.connectance = connectance(web)
        # Specificity and all
        self.generality = generality(web)
        self.vulnerability = vulnerability(web)
        self.specificity = specificity(web)
        self.rr = RR(web)
        self.ssi = ssi(web)
        self.bperf = web.max(1)
        # Nestedness
        self.nodf_strict = nodf_strict
        NODF = nodf(web,strict=nodf_strict,use_c=self.use_c)
        self.nodf = NODF[0]
        self.nodf_up = NODF[2]
        self.nodf_low = NODF[1]
        # Placeholder for modularity
        self.modules = modules(self)
        # Placeholder for tests
        self.tests = test(self)
        # Placeholder for species names
        self.upnames = range(self.upsp)
        self.lonames = range(self.losp)
        # For the name
        self.name = ''
        self.G = nx.DiGraph()
    def __str__(self):
        s = 'Network '+self.name
        s+= '\t['+str(self.losp)+'x'+str(self.upsp)+'] Co = '+str(self.connectance)+'\n'
        return s
    def nxExport(self):
        """
        Expoort the bipartite graph to networkX
        """
        for up_sp in xrange(self.upsp):
            for lo_sp in xrange(self.losp):
                if self.web[up_sp][lo_sp] > 0:
                    t_s = 'top_'+str(self.upnames[up_sp])
                    b_s = 'bot_'+str(self.lonames[lo_sp])
                    l_s = self.web[up_sp][lo_sp]
                    self.G.add_edge(t_s, b_s,weight=l_s)
    def save(self):
        if self.name == '':
            self.name = 'web'
        fname = self.name + '.bipartite'
        file_bip = open(fname, 'w')
        pickle.dump(self, file_bip)
        print 'The bipartite network was saved to '+fname
        return fname
    def txt(self):
        """
        Prints a text version of the network to the console
        Can be used to view quickly the shape of a network
        """
        for line in self.adjacency:
            s = ''
            for char in line:
                if char > 0:
                    s += u'\u2588'
                else:
                    s += '-'
            print s
        return 0
    def plot(self,asnest=True,asbeads=False,colors=True):
        filename = self.name
        if asnest:
            filename += '_nested'
            tW = sortbydegree(self)
            W = bipartite(tW[0])
            W.upnames=tW[1]
            W.lonames=tW[2]
            if asbeads:
                plotBMatrix(W,filename=filename,withcolors=colors)
            else:
                filename += '_matrix'
                plotMatrix(W,filename=filename,withcolors=colors)
        else:
            if not self.modules.done:
                self.modules.detect(self.q_c)
            filename += '_modular'
            g = np.copy(self.modules.up_modules)
            h = np.copy(self.modules.low_modules)
            W = sortbymodule(self,g,h)
            W.modules.N = self.modules.N
            W.modules.up_modules = self.modules.up_modules
            W.modules.low_modules = self.modules.low_modules
            W.modules.Q = self.modules.Q
            W.modules.Qr = self.modules.Qr
            if asbeads:
                plotBModules(W,filename=filename,withcolors=colors)
            else:
                filename += '_matrix'
                plotModules(W,filename=filename,withcolors=colors)
    def networklevel(self,toScreen=True,toFile=True):
        """
        Print the network level statistics to a file / to screen
        """
        fname = 'networklevel.txt'
        fileExt = os.path.exists(fname)
        if fileExt:
            f = open(fname,'a')
        else:
            f = open(fname,'w')
        header = 'name\tCo\tS\tr_up\tr_lo\tnodf\tnodf_up\tnodf_low\taspe\tarr\tresp\tinc'
        info = str(self.name)+'\t'+str(self.connectance)+'\t'+str(self.size)+'\t'+str(self.upsp)+'\t'+str(self.losp)+'\t'+str(self.nodf)+'\t'+str(self.nodf_up)+'\t'+str(self.nodf_low)+ '\t' + str(
            np.mean(self.specificity)) + '\t' + str(np.mean(self.rr)) + '\t' + str(IR(self)[0]) + '\t' + str(IR(self)[1])
        if self.modules.done:
            header += '\tmN\tmQb\tmQr'
            info += '\t'+str(self.modules.N)+'\t'+str(self.modules.Q)+'\t'+str(self.modules.Qr)
        if len(self.tests.devnest) > 0:
            header += '\tnodf_sim\tnodf_pval\tnodf_icLow\tnodf_icUp'
            info += "\t"+str(round(self.tests.devnest[2],2))+"\t"+str(self.tests.devnest[1])+"\t"+str(round(self.tests.devnest[3],2))+"\t"+str(round(self.tests.devnest[4],2))
        if len(self.tests.devnest_lo) > 0:
            header += '\tnodf_low_sim\tnodf_low_pval\tnodf_low_icLow\tnodf_low_icUp'
            info += '\t'+str(round(self.tests.devnest_lo[2],2))+"\t"+str(self.tests.devnest_lo[1])+"\t"+str(round(self.tests.devnest_lo[3],2))+"\t"+str(round(self.tests.devnest_lo[4],2))
        if len(self.tests.devnest_up) > 0:
            header += '\tnodf_up_sim\tnodf_up_pval\tnodf_up_icLow\tnodf_up_icUp'
            info += "\t"+str(round(self.tests.devnest_up[2],2))+"\t"+str(self.tests.devnest_up[1])+"\t"+str(round(self.tests.devnest_up[3],2))+"\t"+str(round(self.tests.devnest_up[4],2))
        if len(self.tests.devqr) > 0:
            header += '\tmQr_sim\tmQr_pval\tmQr_icLow\tmQr_icUp'
            info += "\t"+str(round(self.tests.devqr[2],2))+"\t"+str(self.tests.devqr[1])+"\t"+str(round(self.tests.devqr[3],2))+"\t"+str(round(self.tests.devqr[4],2))
        if len(self.tests.devqb) > 0:
            header += '\tmQb_sim\tmQb_pval\tmQb_icLow\tmQb_icUp'
            info += "\t"+str(round(self.tests.devqb[2],2))+"\t"+str(self.tests.devqb[1])+"\t"+str(round(self.tests.devqb[3],2))+"\t"+str(round(self.tests.devqb[4],2))
        header += '\ttReps\ttMod'
        info += '\t'+str(self.tests.replicates)+'\t'+str(self.tests.model.__name__)
        if toScreen:
            print header
            print info
        if toFile:
            if not fileExt:
                f.write(header)
            f.write('\n')
            f.write(info)
        f.close()
        return 0
    def specieslevel(self,toScreen=True,toFile=True):
        """
        Write the species level informations
        """
        if toFile:
            f = open(self.name+'-sp.txt','w')
        Header = 'sp\tlev\tdeg\tspe\tssi\trr'
        if self.modules.done:
            Header+= '\t'
            Header+= 'MOD'
        if toScreen:
            print Header
        if toFile:
            f.write(Header+'\n')
        for tls in xrange(self.upsp):
            SpInfo = str(self.upnames[tls])+'\t'
            SpInfo+= 'top'+'\t'
            SpInfo+= str(self.generality[tls])+'\t'
            SpInfo+= str(self.specificity[tls])+'\t'
            SpInfo+= str(self.ssi[tls])+'\t'
            SpInfo+= str(self.rr[tls])
            if self.modules.done:
                SpInfo+= '\t'
                SpInfo+= str(self.modules.up_modules[tls])
            if toScreen:
                print SpInfo
            if toFile:
                f.write(SpInfo+'\n')
        for lls in xrange(self.losp):
            SpInfo = str(self.lonames[lls])+'\t'
            SpInfo+= 'bottom'+'\t'
            SpInfo+= str(self.vulnerability[lls])+'\t'
            SpInfo+= '-----'+'\t'
            SpInfo+= '-----'+'\t'
            SpInfo+= '-----'
            if self.modules.done:
                SpInfo+= '\t'
                SpInfo+= str(self.modules.low_modules[lls])
            if toScreen:
                print SpInfo
            if toFile:
                f.write(SpInfo+'\n')
        if toFile:
            f.close()
            print 'Species-level infos for the dataset '+self.name+' were written to '+self.name+'-sp.txt\n'
        return 0

## Sort by modules
def sortbymodule(W,g,h):
    sg = sorted(np.copy(g),reverse=True)
    sh = sorted(np.copy(h),reverse=True)
    # Step 1 : Sort a matrix by module
    ## VOID VECTORS FOR THE SORTED SPECIES NAMES
    oTnames = W.upnames
    oBnames = W.lonames
    vTnames = np.copy(oTnames)
    vBnames = np.copy(oBnames)
    ## Step 1a : sort TLO
    rG = rank(g)
    nW = np.zeros((W.upsp,W.losp))
    for ro in range(0,W.upsp):
        nW[rG[ro]] = W.web[ro]
        vTnames[rG[ro]] = oTnames[ro]
    ## Step 1b : sort BLO
    nW = nW.T
    dW = np.zeros((W.upsp,W.losp)).T
    rG = rank(h)
    for ro in range(0,W.losp):
        dW[rG[ro]] = nW[ro]
        vBnames[rG[ro]] = oBnames[ro]
    web = np.copy(dW.T)
    # New temp files for the names
    oBnames = np.copy(vBnames)
    oTnames = np.copy(vTnames)
    # Step 2 : Sort each module by degree
    uniqueMod = sorted(uniquify(sg),reverse=True)
    ## Step 2a : sort TLO
    totalMadeInt = 0
    tempIntCnt = 0
    tdeg = generality(web)
    nweb = np.zeros(np.shape(web))
    for module in uniqueMod:
        totalMadeInt += tempIntCnt
        tempIntCnt = 0
        cdeg = []
        for sp in range(len(tdeg)):
            if sg[sp] == module:
                cdeg.append(tdeg[sp])
        rnk = rank(cdeg)
        for ro in range(len(rnk)):
            nweb[totalMadeInt+rnk[ro]] = web[totalMadeInt+ro]
            vTnames[totalMadeInt+rnk[ro]] = oTnames[totalMadeInt+ro]
            tempIntCnt += 1
    web = np.copy(nweb.T)
    ## Step 2b : sort BLO
    totalMadeInt = 0
    tempIntCnt = 0
    tdeg = generality(web)
    nweb = np.zeros(np.shape(web))
    for module in uniqueMod:
        totalMadeInt += tempIntCnt
        tempIntCnt = 0
        cdeg = []
        for sp in range(len(tdeg)):
            if sh[sp] == module:
                cdeg.append(tdeg[sp])
        rnk = rank(cdeg)
        for ro in range(len(rnk)):
            nweb[totalMadeInt+rnk[ro]] = web[totalMadeInt+ro]
            vBnames[totalMadeInt+rnk[ro]] = oBnames[totalMadeInt+ro]
            tempIntCnt += 1
    Fweb = bipartite(np.copy(nweb.T))
    Fweb.upnames = vTnames
    Fweb.lonames = vBnames
    return Fweb
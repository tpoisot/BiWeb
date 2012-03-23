import pickle
import os
import tempfile
import tkFileDialog

from ..bipartite_class import *

def load(fname):
    file_bip = open(fname, 'r')
    bip = pickle.load(file_bip)
    file_bip.close()
    return bip

def openWeb(file='',t=False,name='',species_names=False):
    if species_names:
        web = oNW(file,t,name)
    else:
        web = oUW(file,t,name)
    return web

def oUW(file='',t=False,name=''):
    if file == '':
        filename = tkFileDialog.askopenfilename()
    else:
        filename = file
        # Read the web
    w = bipartite(readweb(filename),t=t)
    w.name = name
    return w


def oNW(file='',t=False,name=''):
    if file == '':
        filename = tkFileDialog.askopenfilename()
    else:
        filename = file
    upnames = []
    lonames = []
    f = tempfile.NamedTemporaryFile(delete=False)
    # Read the web
    fweb = open(filename,'r')
    cline = 0
    for line in fweb:
        if cline == 0:
            lonames = line.split()
        else:
            spl_line = line.split()
            upnames.append(spl_line[0])
            tintmat= []
            for i in range(1,len(spl_line)):
                f.write(str(float(spl_line[i]))+' ')
            f.write('\n')
        cline += 1
    f.close()
    web = oUW(f.name,t,name)
    os.unlink(f.name)
    web.upnames=upnames
    web.lonames=lonames
    return web

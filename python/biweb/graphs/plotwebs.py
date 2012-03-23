import scipy as sp
import numpy as np

from ..bipartite_class import *
from pyx import *
from ..base import *
from ..mod import *

# Plot a modular web as beads
def plotBModules(w,filename='web',withcolors=True):
    GS = 0.5
    # Organise web by communities
    W = np.copy(w.web)
    g = np.copy(w.modules.up_modules)
    cg = sorted(g,reverse=True)
    h = np.copy(w.modules.low_modules)
    ch = sorted(h,reverse=True)
    # Aspect ratio (these parameters seem to be OK)
    yp = max((max(w.upsp,w.losp)/float(min(w.upsp,w.losp))*5),8)
    # Next...
    ListOfColors = [color.gradient.Hue.select(i, (w.modules.N+1)) for i in range(w.modules.N+1)]
    ###### Plot
    # Get the line length
    MLink = max(w.bperf)
    # Get the positions of the nodes
    mnodes = max(w.upsp,w.losp)
    top_height = spread(range(w.upsp),1,mnodes)
    bot_height = spread(range(w.losp),1,mnodes)
    # Plot
    c = canvas.canvas()
    # Plot the lines first
    for i in range(w.upsp):
        for j in range(w.losp):
            if W[i][j] > 0:
                lwid = round((W[i][j]/float(MLink)),0)*0.05
                if cg[i] == ch[j]:
                    continue
                else:
                    c.stroke(path.line(1, top_height[i], yp, bot_height[j]),[deco.stroked([color.cmyk.Gray]),style.linewidth(lwid)])

    # Plot the lines first
    for i in range(w.upsp):
        for j in range(w.losp):
            if W[i][j] > 0:
                lwid = round((W[i][j]/float(MLink)),0)*0.05
                if cg[i] == ch[j]:
                    c.stroke(path.line(1, top_height[i], yp, bot_height[j]),[deco.stroked([color.gray.black]),style.linewidth(lwid)])
                else:
                    continue
    # Plot the beads
    for i in range(w.upsp):
        c.text(0.5, top_height[i], w.upnames[i],[text.parbox(2), text.halign.right])
        if withcolors:
            CCol = ListOfColors[(cg[i]-1)]
        else:
            CCol = color.gray.black
        c.fill(path.circle(1, top_height[i], GS*0.6),[deco.stroked.clear,deco.filled([CCol])])
    for j in range(w.losp):
        c.text(yp+0.5, bot_height[j], w.lonames[j],[text.parbox(2), text.halign.left])
        if withcolors:
            CCol = ListOfColors[(ch[j]-1)]
        else:
            CCol = color.gray.black
        c.fill(path.circle(yp, bot_height[j], GS*0.6),[deco.stroked.clear,deco.filled([CCol])])
    c.writePDFfile(filename)
    return 0


# Plot a nested web as beads
def plotBMatrix(w,filename='web',withcolors=True):
    GS = 0.5
    W = np.copy(w.web)
#	w = bipartite(w)
    # Aspect ratio (these parameters seem to be OK)
    yp = max((max(w.upsp,w.losp)/float(min(w.upsp,w.losp))*5),8)
    ###### Plot
    # Get the line length
    MLink = max(w.bperf)
    # Get the positions of the nodes
    mnodes = max(w.upsp,w.losp)
    top_height = spread(range(w.upsp),1,mnodes)
    bot_height = spread(range(w.losp),1,mnodes)
    # Plot
    c = canvas.canvas()
    # Plot the lines first
    for i in range(w.upsp):
        for j in range(w.losp):
            if W[i][j] > 0:
                lwid = round((W[i][j]/float(MLink)),0)*0.05
                c.stroke(path.line(1, top_height[i], yp, bot_height[j]),[deco.stroked([color.gray.black]),style.linewidth(lwid)])
    # Plot the beads
    for i in range(w.upsp):
        c.text(0.5, top_height[i], w.upnames[i],[text.parbox(2), text.halign.right])
        c.fill(path.circle(1, top_height[i], GS*0.6),[deco.stroked.clear,deco.filled([color.gray.black])])
    for j in range(w.losp):
        c.text(yp+0.5, bot_height[j], w.lonames[j],[text.parbox(2), text.halign.left])
        c.fill(path.circle(yp, bot_height[j], GS*0.6),[deco.stroked.clear,deco.filled([color.gray.black])])
    c.writePDFfile(filename)
    return 0


# Plot a web as a nested matrix
def plotMatrix(w,filename='web',withcolors=True):
    GS = 0.5
    W = np.copy(w.web)
#	w = bipartite(w)
    c = canvas.canvas()
    if withcolors:
        MLink = max(w.bperf)
        for i in range(w.upsp):
            for j in range(w.losp):
                W[i][j] = round((W[i][j]/float(MLink))*100,0)
        # Define the color gradient
        ListOfColors = [color.gradient.Gray.select(i, 101) for i in range(101)]
    for i in range(w.upsp):
        c.text(-0.1, GS*(i+0.7), w.upnames[i],[text.parbox(2), text.halign.right])
        for j in range(w.losp):
            if i == 0:
                c.text(GS*(j+0.9), -0.1, w.lonames[j],[text.parbox(2), text.halign.right, trafo.rotate(90)])
            if W[i][j] > 0:
                xc = GS*j+(GS/2)
                yc = GS*i+(GS/2)
                xd = 0.8*GS
                yd = 0.8*GS
                if withcolors:
                    c.stroke(path.rect(xc, yc, xd, yd),[deco.stroked.clear,deco.filled([ListOfColors[int(W[i][j])]])])
                else:
                    c.stroke(path.rect(xc, yc, xd, yd),[deco.stroked.clear,deco.filled([color.gray.black])])
    c.writePDFfile(filename)
    return 0


# Plot a web as a modular matrix
def plotModules(w,filename='web',withcolors=True):
    GS = 0.5
    # Organise web by communities
    W = np.copy(w.web)
    g = np.copy(w.modules.up_modules)
    cg = sorted(g,reverse=True)
    h = np.copy(w.modules.low_modules)
    ch = sorted(h,reverse=True)
#	w = mini_bipartite(w)
    ListOfColors = [color.gradient.Hue.select(i, (w.modules.N+1)) for i in range(w.modules.N+1)]
    # Plot
    c = canvas.canvas()
    for i in range(w.upsp):
        c.text(-0.1, GS*(i+0.7), w.upnames[i],[text.parbox(2), text.halign.right])
        for j in range(w.losp):
            if i == 0:
                c.text(GS*(j+0.9), -0.1, w.lonames[j],[text.parbox(2), text.halign.right, trafo.rotate(90)])
            if W[i][j] > 0:
                xc = GS*j+(GS/2)
                yc = GS*i+(GS/2)
                xd = 0.8*GS
                yd = 0.8*GS
                if cg[i] == ch[j]:
                    if withcolors:
                        CCol = ListOfColors[(cg[i]-1)]
                    else:
                        CCol = color.gray.black
                else:
                    CCol = color.cmyk.Gray
                c.stroke(path.rect(xc, yc, xd, yd),[deco.stroked.clear,deco.filled([CCol])])
    c.writePDFfile(filename)
    return 0

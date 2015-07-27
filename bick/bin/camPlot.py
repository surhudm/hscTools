#!/usr/bin/env python
# Original filename: camPlot.py
#
# Author: Steve Bickerton
# Email: 
# Date: Tue 2013-12-31 13:02:29
# 
# Summary: 
# 

import sys
import os
import re
import math
import argparse
import numpy
import datetime
import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas

import lsst.pex.policy           as pexPolicy
import lsst.afw.cameraGeom       as afwCG
import lsst.afw.cameraGeom.utils as afwCGU

import hsc.tools.bick.utils           as hscUtil


#############################################################
#
# Main body of code
#
#############################################################

def main(infile, ccds=None, camera="hsc", cmap="copper", cols=[0,1],
         nsig=3.0, percent=False, textcolor='k', out=None):

    ###########################
    # build the camera
    policy_file = {
        'hsc': os.path.join(os.getenv("OBS_SUBARU_DIR"), "hsc", "hsc_geom.paf"),
        'sc' : os.path.join(os.getenv("OBS_SUBARU_DIR"), "suprimecam", "Full_Suprimecam_geom.paf")
        }
    geomPolicy = afwCGU.getGeomPolicy(policy_file[camera.lower()])
    camera     = afwCGU.makeCamera(geomPolicy)

    
    ###########################
    # load the data
    data = {}
    if infile:
        load = numpy.loadtxt(infile)
        vals = load[:,cols[1]]
        for i in range(len(vals)):
            ccd = int(load[i,cols[0]])
            if len(cols) == 3:
                amp = int(load[i,cols[2]])
            else:
                amp = 0
            if ccd not in data:
                data[ccd] = {}
            data[ccd][amp] = vals[i]
    else:
        vals = []
        for r in camera:
            for c in afwCG.cast_Raft(r):
                ccd = afwCG.cast_Ccd(c)
                ccdId = ccd.getId().getSerial()
                data[ccdId] = {}
                val = 1.0
                if ccdId > 103:
                    val = 0.8
                for a in ccd:
                    amp = afwCG.cast_Amp(a)
                    ampId = amp.getId().getSerial() - 1
                    data[ccdId][ampId] = val
                    vals.append(val)
                if len(data[ccdId]) == 0:
                    data[ccdId][0] = val
                    vals.append(val)
        vals = numpy.array(vals)
        
    mean = vals.mean()
    med  = numpy.median(vals)
    std  = vals.std()
    vmin, vmax = med - nsig*std, med+nsig*std
    

    ###########################
    # make the plot
    fig = figure.Figure(figsize=(7,7))
    canvas = FigCanvas(fig)

    if infile:
        rect = (0.06, 0.12, 0.76, 0.76)
    else:
        rect = (0.06, 0.06, 0.88, 0.88)
        
    fpa_fig = hscUtil.FpaFigure(fig, camera, rect=rect)
    
    for i_ccd, amplist in data.items():

        hide = False
        if ccds and (i_ccd not in ccds):
            hide = True
        
        ax = fpa_fig.getAxes(i_ccd)
        fpa_fig.highlightAmp(i_ccd, 0)
        nq = fpa_fig.detectors[i_ccd].getOrientation().getNQuarter()
        nx, ny = 4, 8
        if nq % 2:
            ny, nx = nx, ny
        firstAmp = sorted(amplist.keys())[0]
        im = numpy.zeros((ny, nx)) + amplist[firstAmp]

        print i_ccd
        for amp, val in amplist.items():
            useVal = val
            if hide:
                useVal = 0.0
            if nq == 0:
                im[:,amp] = useVal
            if nq == 1 or nq == -3:
                im[3-amp,:] = useVal
            if nq == 2:
                im[:,3-amp] = useVal
            if nq == -1 or nq == 3:
                im[amp,:] = useVal
            
        im_ax = ax.imshow(im, cmap=cmap, vmax=vmax, vmin=vmin, interpolation='nearest')

        fpa_fig.addLabel(i_ccd, [str(i_ccd)], color=textcolor)


    #############################
    # the colorbar
    ylo, yhi = vmin, vmax

    if infile:
        rect = (0.91, 0.25, 0.02, 0.6)
        
        # real units
        cax  = fig.add_axes(rect)        
        cax.get_yaxis().get_major_formatter().set_useOffset(False)        
        cax.set_ylim([ylo, yhi])
        cax.get_xaxis().set_ticks([])
        cbar = fig.colorbar(im_ax, cax=cax)

        # mirror the values.  std or in percent (fractional) if requested
        caxp = cax.twinx()
        caxp.get_xaxis().set_ticks([])
        caxp.get_yaxis().get_major_formatter().set_useOffset(False)        
        if percent:
            caxp.set_ylim([(med-ylo)/(yhi-ylo), (yhi-med)/(yhi-ylo)])
        else:
            caxp.set_ylim([-nsig*std, nsig*std])


        for t in cax.get_yticklabels():
            t.set_size('small')
        for t in caxp.get_yticklabels():
            t.set_size("small")


    #################################
    # add stats
    if infile:
        stats = {"Mean":mean, "Med":med, "Std":std, "Max":vals.max(), "Min":vals.min()}
        order = ["Mean", "Med", "Std", "Min", "Max"]
        i = 0
        for k in order:
            v = stats[k]
            ax.text(0.8, 0.03+0.02*i, "%-10s %.3g"%(k+":",v), fontsize=9,
                    horizontalalignment='left', verticalalignment='center', transform=fig.transFigure)
            i += 1

    #################################
    # title and write it
    date = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    title = "File: %s  Cols: %s  %s" % (infile, ":".join([str(x+1) for x in cols]), date) if infile else ""
    fig.suptitle(title)

    outfile = out if out else "camview-%s.png" % infile
    fig.savefig(outfile)


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action='version', version="0.0")
    parser.add_argument("infile", type=str, help="")
    parser.add_argument("-C", "--camera", type=str, default="hsc", choices=("hsc","sc","lsst"),  help="")
    parser.add_argument("-m", "--cmap", type=str, default="copper",
                        choices=("copper","gray","jet", "prism"),  help="")
    parser.add_argument("-d", "--datacols", type=str, default="1:2", help="Columns format='CCD:datavalue'")
    parser.add_argument("-n", "--nsig", type=float, default=3.0, help="N standard devs for cmap")
    parser.add_argument("-p", "--percent", default=False, action='store_true',
                        help="display percent on complementary colorbar")
    parser.add_argument("-c", "--ccds", default=None, help="Ccds to show")
    parser.add_argument("-t", "--textcolor", default='k', help="Color for labels")
    parser.add_argument("-o", "--out", default=None, help="Output plot name")
    args = parser.parse_args()

    infile = None if args.infile == 'none' else args.infile
    
    cols = [int(x)-1 for x in args.datacols.split(":")]
    ccds=None
    if args.ccds:
        ccds = set([int(x) for x in hscUtil.idSplit(args.ccds)])
    
    main(infile, ccds=ccds, camera=args.camera, cmap=args.cmap, cols=cols,
         nsig=args.nsig, percent=args.percent, textcolor=args.textcolor, out=args.out)

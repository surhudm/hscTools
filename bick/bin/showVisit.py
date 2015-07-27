#!/usr/bin/env python
# Original filename: showVisit.py
#
# Author: Steve Bickerton
# Email: 
# Date: Thu 2013-12-12 15:01:41
# 
# Summary: 
# 

import sys
import os
import re
import math
import argparse
import datetime
import lsst.daf.persistence          as dafPersist
import hsc.pipe.base.butler          as hscButler
import numpy
import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
import matplotlib.font_manager as fm

import hsc.tools.bick.utils as hscUtil

import lsst.afw.cameraGeom          as camGeom
import lsst.afw.cameraGeom.utils    as camGeomUtils


scaleLookup = {
    'linear' : lambda x : x,
    'histeq' : hscUtil.histeq,
    }

def compute(im1, im2, op):
    if op == 'p':
        ret = im1 + im2
    if op == 'm':
        ret = im1 - im2
    if op == 'd':
        ret = im1 / im2
    if op == 't':
        ret = im1 * im2
    return ret


calibTypes = "dark", "bias", "flat"

#############################################################
#
# Main body of code
#
#############################################################

def main(rerun, visit1, visit2, op, ccdno,
         datatype='calexp', scale=None, root=None, invert=False, cmap='gray',
         vmax=None, annotate=None, bins=16, rerun2=None, showCbar=False, vsig=None,
         showAnnotate=False, percent=None, hilite=None):

    vsigDefault = 5.0
    visit1 = int(visit1)
    if visit2:
        visit2 = int(visit2)
    ccdno = set([int(x) for x in hscUtil.idSplit(ccdno)])
    
    butler1 = hscUtil.getButler(rerun, root=root)
    if rerun2:
        butler2 = hscUtil.getButler(rerun2, root=root)
    else:
        butler2 = butler1

    if datatype not in calibTypes:
        dataIds1 = butler1.queryMetadata(datatype, "ccd", format=["visit", "ccd"], dataId={'visit':visit1})
        dataIds1 = [{'visit':x[0], 'ccd':x[1]} for x in dataIds1]
    else:
        dataIds1 = [{'visit': visit1, 'ccd':x} for x in range(104)]
        
    dataRef1 = None
    dataRef2 = None
    
    if visit2:
        dataIds2 = butler2.queryMetadata(datatype, "ccd", format=["visit", "ccd"], dataId={'visit':visit2})
        dataIds2 = [{'visit':x[0], 'ccd':x[1]} for x in dataIds2]
    
    # flip the color map
    if invert:
        cmap = re.sub("_r$", "", cmap) if re.match('_r$', cmap) else '_r'

    # sleezy, but if vmax isn't set and we're linear, just use ccd in the middle to norm 
    if vmax is None and scale != 'histeq':
        vmax = 'c049'

    # handle the gray scale normalization
    vmin = 0.0
    mdRef1 = None
    mdRef2 = None
    if vmax:
        if scale == 'histeq':
            raise ValueError("Cannot specify vmax with histeq scaling.")
        
        # if it identifies a CCD to use
        if re.match("^c", vmax):
            vmaxCcd = int(re.sub("c", "", vmax))
            try:
                dataRef1 = hscButler.getDataRef(butler1, {'visit':visit1, 'ccd':vmaxCcd})
                imgRef1  = dataRef1.get(datatype).getMaskedImage().getImage().getArray()
                if datatype not in calibTypes:
                    mdRef1   = dataRef1.get(datatype+'_md', immediate=True)
                if visit2:
                    dataRef2 = hscButler.getDataRef(butler2, {'visit':visit2, 'ccd':vmaxCcd})
                    imgRef2  = dataRef2.get(datatype).getMaskedImage().getImage().getArray()
                    mdRef2   = dataRef2.get(datatype+'_md', immediate=True)
                    img_op  = hscUtil.rebin(compute(imgRef1, imgRef2, op), bins)
                    med     = numpy.median(img_op)
                    std     = numpy.std(img_op)
                    if not vsig:
                        vsig = vsigDefault
                    delta = vsig*std
                    if percent:
                        delta = percent*med
                    vmin    = med - abs(delta)
                    vmax    = med + abs(delta)
                else:
                    if showAnnotate:
                        exp1 = dataRef1.get(datatype)
                        aval   = float(exp1.getMetadata().get(annotate))
                        med = aval
                        delta = vsig
                        if not vsig:
                            delta = 0.5*aval
                        if percent:
                            delta = percent*aval
                        vmin    = med - abs(delta)
                        vmax    = med + abs(delta)
                    else:
                        img_op  = imgRef1
                        med     = numpy.median(img_op)
                        sig     = numpy.sqrt(med)
                        delta   = vsigDefault*sig
                        if vsig:
                            delta = vsig*sig
                        if percent:
                            delta = percent*med
                        vmin    = med - abs(delta)
                        vmax    = med + abs(delta)
                        if not vsig and not percent:
                            vmin    = 0.5*med
                            vmax    = med + 5.0*sig

                        
            except Exception, e:
                raise RuntimeError("Could not get stats on vmax CCD" + str(vmax)+ "  Exiting." + str(e))

        elif re.search(":", vmax):
            vmaxCcd = None
            vmin, vmax = [float(x) for x in vmax.split(":")]
            med = 0.5*(vmax + vmin)
        else:
            vmaxCcd = None
            vmax = float(vmax)
            med = 0.5*(vmax + vmin)
            
            
    ###########################
    fig = figure.Figure((8,8))
    canvas = FigCanvas(fig)

    l, b, w, h = 0.08, 0.12, 0.84, 0.78
    rect = l, b, w, h
    if showCbar:
        rect = (0.06, 0.12, 0.76, 0.76)
    fpa_fig = hscUtil.FpaFigure(fig, butler1.get('camera'), rect=rect)

    im_ax = None
    i_show = 0
    for dataId1 in dataIds1:
        if dataId1['ccd'] not in ccdno:
            continue
        print dataId1,
        
        if visit2:
            dataId2 = {'visit':visit2, 'ccd':dataId1['ccd']}
            print dataId2
        else:
            print ""
            
        try:
            dataRef1 = hscButler.getDataRef(butler1, dataId1)
            exp1 = dataRef1.get(datatype)
            img1 = None
            if not showAnnotate:
                img1 = exp1.getMaskedImage().getImage().getArray()
                img1 = img1
            if visit2:
                dataRef2 = hscButler.getDataRef(butler2, dataId2)
                exp2 = dataRef2.get(datatype)
                img2 = exp2.getMaskedImage().getImage().getArray()
        except Exception, e:
            #raise
            print "getDataRef() failed for ", visit1, visit2, ccdno, str(e)
            continue

        ax = fpa_fig.getAxes(dataId1['ccd'])


        # labels
        labels = [ str(dataId1['ccd']) ]
        aValues = None
        if annotate:
            for a in annotate.split('|'):
                if a == 'NQUARTER':
                    labels.append(str(exp.getDetector().getOrientation().getNQuarter()))
                if a in exp1.getMetadata().paramNames():
                    aval = exp1.getMetadata().get(annotate)
                    labels.append(aval)
                    if not aValues:
                        aValues = aval
        fpa_fig.addLabel(dataId1['ccd'], labels)

        if hilite is not None:
            fpa_fig.highlightAmp(dataId1['ccd'], hilite)
            
            
        imshow_kwargs = {}
        if scale == 'linear':
            imshow_kwargs = {'vmin': vmin, 'vmax': vmax}
        if showAnnotate:
            ny, nx = 4, 2
            if exp1.getDetector().getOrientation().getNQuarter() % 2:
                ny, nx = 2, 4
            imtmp = numpy.ones((ny, nx))*float(aval)
        else:
            if visit2:
                img_op = compute(img1, img2, op)
            else:
                img_op = img1
            # scale as requested
            img_op = scaleLookup[scale](img_op)
            imtmp = hscUtil.rebin(img_op, bins)[::-1]
        
        im_ax = ax.imshow(imtmp, cmap=cmap, **imshow_kwargs)

        

        
    if showCbar and im_ax and scale != "histeq":
        ylo, yhi = vmin, vmax
        rect = (0.91, 0.2, 0.02, 0.6)

        cax  = fig.add_axes(rect)        
        cax.set_ylim([ylo, yhi])
        cax.get_xaxis().set_ticks([])
        cax.get_yaxis().get_major_formatter().set_useOffset(False)        
        cbar = fig.colorbar(im_ax, cax=cax)

        caxp = cax.twinx()
        caxp.set_ylim([ylo/med, yhi/med])
        caxp.get_xaxis().set_ticks([])
        caxp.get_yaxis().get_major_formatter().set_useOffset(False)        
        
        for t in cax.get_yticklabels():
            t.set_size('small')
        for t in caxp.get_yticklabels():
            t.set_size("small")

        
    i = 0
    # get some header info to print labels
    hdrcards = ["OBJECT", "FILTER", "EXPTIME", "DATE-OBS"]
    #ax = fig.add_axes((l, 0.05, w, 0.08))
    # kill the ticks
    #for tick in ax.get_xticklines() + ax.get_yticklines() + ax.get_yticklabels() + ax.get_xticklabels():
    #    tick.set_visible(False)
    ops = {'p': 'plus', "m": 'minus', 'd': 'div', "t": 'times'}
    for card in hdrcards:
        if not mdRef1:
            if dataRef1 and datatype not in calibTypes:
                mdRef1 = dataRef1.get(datatype+'_md', immediate=True)
        if visit2 and not mdRef2 and datatype not in calibTypes:
            mdRef2 = dataRef2.get(datatype+'_md', immediate=True)
            
        if mdRef1 and card in mdRef1.paramNames():
            ax.text(0.25, 0.02+0.015*i, card+": "+str(mdRef1.get(card)), fontsize=8,
                    horizontalalignment='left', verticalalignment='center', transform=fig.transFigure)
            if visit2:
                ax.text(0.5, 0.04, ops[op], fontsize=8,
                        horizontalalignment='center', verticalalignment='center', transform=fig.transFigure)
                ax.text(0.6, 0.02+0.015*i, card+": "+str(mdRef2.get(card)), fontsize=8,
                        horizontalalignment='left', verticalalignment='center', transform=fig.transFigure)
            i += 1

    rerunStr = re.sub("/", "_", rerun)
    date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if visit2:
        if rerun2:
            rerun2Str = re.sub("/", "_", rerun2)
            key = "%s_%s-%d%s%d-%s" % (rerunStr, rerun2Str, visit1, op, visit2, scale[0].upper())
        else:
            key = "%s-%d%s%d-%s" % (rerunStr, visit1, op, visit2, scale[0].upper())
    else:
        key = "%s-%d-%s" % (rerunStr, visit1, scale[0].upper())
    if showAnnotate:
        key += "-%s" % (annotate.split("|")[0])
    fig.suptitle(key + " (" + datatype + ") "+ date)
    fig.savefig("fpa-%s.png" % (key))

    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action='version', version="0.0")
    parser.add_argument("rerun", type=str, help="")
    parser.add_argument("visit", type=str, help="")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--plus", default=None)
    group.add_argument("--minus", default=None)
    group.add_argument("--div", default=None)
    group.add_argument("--times", default=None)
    
    parser.add_argument("ccd", type=str, help="")
    parser.add_argument("-a", "--annotate", default=None, help="Header values to display.  Options may include: T_CCDTV, T_GAIN1, GAINEFF, SKYLEVEL, SKYSIGMA, HIERARCH FLATNESS_PP, DET-TMED, HIERARCH NOBJ_BRIGHT,HIERARCH runEndCpuTime, HIERARCH fwhmRobust")
    parser.add_argument("-A", "--showAnnotate", default=False, action='store_true', help="Use the value of the first annotate term to color the figure.  Use '-V' to control the cmap range above/below the ")
    parser.add_argument("-b", '--bins', type=int, default=16, help="Binning to use.")
    parser.add_argument("-C", "--nocbar", default=False, action='store_true', help="Do *NOT* show color bar")
    parser.add_argument("-c", "--cmap", type=str, default='copper',
                        choices=("gray", 'jet', 'copper', "gist_rainbow", "Spectral", "Accent", "Blues",
                                 "BrBG", "BuPu", "Dark2", "GnBu", "Greens", "Greys", "OrRd", "Oranges",
                                 "PRGn", "Paired", "Pastel1", "Pastel2", "PiYG", "PuBu", "PuBuGn", "PuOr",
                                 "PuRd", "Purples", "RdBu", "RdGy", "RdPu", "RdYlBu", "RdYlGn", "Reds",
                                 "Set1", "Set2", "Set3", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd", "autumn",
                                 "binary", "bone", "cool", "flag", "gist_earth", "gist_gray", "gist_heat",
                                 "gist_ncar", "gist_rainbow", "gist_stern", "gist_yarg", "gray", "hot",
                                 "hsv", "pink", "prism", "spectral", "sprint", "summer", "winter"),
                        help="Specify a color map")
    parser.add_argument("-d", "--datatype", default="calexp", help="Type of image to try to load",
                        choices=("raw", "calexp", "postISRCCD", 'dark', 'bias', 'flat'))
    parser.add_argument("-H", '--hilite', type=int, default=None, help="Hilight amp by number.")
    parser.add_argument("-p", "--percent", default=None, type=float,
                        help="Override vsig with a range in percent")
    parser.add_argument("-R", "--rerun2", default=None, help="Rerun for visit2 if different from visit1")
    parser.add_argument("-s", "--scale", type=str, default='linear',  help="Gray scaling.",
                        choices=("histeq", "linear"))
    parser.add_argument("-i", "--invert", default=False, action='store_true', help="invert the gray scale")
    parser.add_argument("-r", "--root", type=str, default=None, help="")
    parser.add_argument("-V", "--vsig", type=float, default=None, help="Sigma for color map")
    parser.add_argument("-x", "--vmax", type=str, default=None, help="Scaling max value.")

    
    args = parser.parse_args()

    visit2 = None
    op = None
    if args.plus:
        visit2 = args.plus
        op = 'p'
    if args.minus:
        visit2 = args.minus
        op = 'm'
    if args.div:
        visit2 = args.div
        op = 'd'
    if args.times:
        visit2 = args.times
        op = 't'
        
    if visit2 and args.showAnnotate:
        print "Cannot use showAnnotate with math operation ... yet."
        sys.exit(1)

    
    main(args.rerun, args.visit, visit2, op, args.ccd, scale=args.scale, root=args.root,
         datatype=args.datatype, invert=args.invert, cmap=args.cmap, vmax=args.vmax,
         annotate=args.annotate, bins=args.bins, rerun2=args.rerun2,
         showCbar=(not args.nocbar), vsig=args.vsig,
         showAnnotate=args.showAnnotate, percent=args.percent, hilite=args.hilite)

    

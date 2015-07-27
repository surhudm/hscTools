#!/usr/bin/env python
# Original filename: visitOverlap.py
#
# Author: Steve Bickerton
# Email: 
# Date: Thu 2013-05-23 12:00:44
# 
# Summary: 
# 

import sys
import os
import re
import argparse
import numpy
import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib.patches              import Rectangle
from matplotlib.patches              import Polygon
from matplotlib.collections          import PatchCollection
from matplotlib                      import cm

import lsst.afw.image                as afwImage
import lsst.afw.coord                as afwCoord
import lsst.afw.geom                 as afwGeom
import lsst.daf.persistence          as dafPersist

import hsc.pipe.base.butler          as hscButler

import hsc.tools.bick.utils as hscUtil

def visitCcdToDataId(visit, ccd):
    visits = hscUtil.idSplit(visit)
    ccds   = hscUtil.idSplit(ccd)
    dataIds = []
    for v in visits:
        for c in ccds:
            dataId = {'visit':v, 'ccd':c}
            dataIds.append(dataId)
    return dataIds


#############################################################
#
# Main body of code
#
#############################################################

def main(visit, rerun, out, ccd="0..103", root=None, 
         point=None, minchips=1, showLabels=False, boresightOnly=False, edge=0.9):

    pcen, rad, pcenVisit = None, None, None
    if point:
        pcenArgs = map(float, point.split(":"))
        if len(pcenArgs) == 2:
            pcenVisit, rad = pcenArgs
        elif len(pcenArgs) == 3:
            cx, cy, rad = pcenArgs
            pcen = afwCoord.Fk5Coord(afwGeom.Point2D(cx, cy))
        else:
            raise ValueError("Unable to parse center 'point', must be ra:dec:rad or visit:rad")

        
    print visit, rerun
    butler = hscUtil.getButler(rerun, root=root)

    visitWidths = []
    visitHeights = []
    
    dataIds = visitCcdToDataId(visit, ccd)
    corners = {}
    boresights = {}
    for dataId in dataIds:
        visit, ccd = dataId['visit'], dataId['ccd']

        if boresightOnly and visit in boresights:
            continue
        
        print "Getting ", visit, ccd
        try:
            dataRef = hscButler.getDataRef(butler, dataId)
        except Exception, e:
            print "getDataRef() failed for ", visit, ccd, str(e)
            continue

        calfile = dataRef.get("calexp_filename", immediate=True)[0]
        noCal = False
        if not os.path.exists(calfile):
            print calfile+" missing."
            noCal = True
        else:
            # get the calexp
            calmd = dataRef.get("calexp_md", immediate=True)
            wcs = afwImage.makeWcs(calmd)
            w, h = calmd.get("NAXIS1"), calmd.get("NAXIS2")

            boresights[dataId['visit']] = [calmd.get("RA2000"), calmd.get("DEC2000")]
            
            # get the corners in RA,Dec
            ll = wcs.pixelToSky(afwGeom.Point2D(1,1)).toFk5()
            lr = wcs.pixelToSky(afwGeom.Point2D(w,1)).toFk5()
            ul = wcs.pixelToSky(afwGeom.Point2D(1,h)).toFk5()
            ur = wcs.pixelToSky(afwGeom.Point2D(w,h)).toFk5()
            
        if noCal:

            rawfile = dataRef.get("raw_filename", immediate=True)[0]
            if not os.path.exists(rawfile):
                print rawfile+" missing too.  continuing."
                continue
            else:
                rawmd = dataRef.get("raw_md", immediate=True)
                raS, decS = rawmd.get("RA2000"), rawmd.get("DEC2000")
                boresights[dataId['visit']] = [raS, decS]
                #ra, dec = c.getRa().asDegrees(), c.getDec().asDegrees()
                crv1, crv2 = map(float, [rawmd.get("CRVAL1"), rawmd.get("CRVAL2")])
                crp1, crp2 = map(float, [rawmd.get("CRPIX1"), rawmd.get("CRPIX2")])

                useDelt = True
                try:
                    cdelt1, cdelt2  = map(float, [rawmd.get("CDELT1"), rawmd.get("CDELT1")])
                except:
                    useDelt= False
                if not useDelt:
                    cd11, cd12, cd21, cd22 = map(float, [rawmd.get("CD1_1"),
                                                         rawmd.get("CD1_2"),
                                                         rawmd.get("CD2_1"),
                                                         rawmd.get("CD2_2")])
                    cdelt1 = 3600*(cd11 + cd12)/206265
                    cdelt2 = 3600*(cd21 + cd22)/206265

                    #cdelt1, cdelt2 = 1.0, 1.0

                nx, ny  = map(int, [rawmd.get("NAXIS1"), rawmd.get("NAXIS2")])

                def trans(x, y):
                    ra = crv1 + (crp1-x)*cdelt1
                    dec = crv2 + (crp2-y)*cdelt2
                    return ra, dec

                ll = afwCoord.Fk5Coord(afwGeom.Point2D(*trans(1, 1)))
                lr = afwCoord.Fk5Coord(afwGeom.Point2D(*trans(nx, 1)))
                ul = afwCoord.Fk5Coord(afwGeom.Point2D(*trans(1, ny)))
                ur = afwCoord.Fk5Coord(afwGeom.Point2D(*trans(nx, ny)))

        if pcenVisit and not pcen:
            if int(pcenVisit) == int(visit):
                cx, cy = boresights[pcenVisit]
                pcen = afwCoord.Fk5Coord(cx, cy)

        if visit not in corners:
            corners[visit] = []
            
        corners[visit].append([ll, lr, ur, ul, ccd])

            
    ccdWidths = []
    ccdHeights = []
    for visit, cornerList in corners.items():

        for bbox in cornerList:
            ll, lr, ul, ur, ccd = bbox

            keep = True
            if pcen:

                # if only doing boresights, don't bother checking the corners
                if boresightOnly:
                    cx, cy = boresights[visit]
                    sep = pcen.angularSeparation(afwCoord.Fk5Coord(cx, cy))
                    if sep.asDegrees() > rad:
                        keep = False
                else:
                    for c in ll,lr,ul,ur:
                        sep = pcen.angularSeparation(c)
                        if sep.asDegrees() > rad:
                            keep = False

            if keep:
                ccdWidths.append(ll.getRa().asDegrees() - lr.getRa().asDegrees())
                ccdHeights.append(ll.getDec().asDegrees() - ul.getDec().asDegrees())
            else:
                if boresights:
                    del boresights[visit]
                del corners[visit]
            
        
    for v, clist in corners.items():
        if pcen and len(clist) < minchips:
            del corners[v]
    visits = sorted(corners.keys())
            
    print "Keeping: "
    for v in visits:
        print " ", v, boresights[v]


    ccdMeanWidth = 0.1
    if ccdWidths:
        ccdMeanWidth = numpy.mean(ccdWidths)
    
    patches = []
    i = 0
    n = len(corners.keys())
    colors = []
    ras = []
    decs = []
    ccd_labels = []
    patch_labels = []
    patch_proxies = []
    proxy_colors = []
    for v in visits:
        L = corners[v]
        color = 1.0
        if n > 1:
            color = 1.0*i/(n-1.0)
        patch_labels.append(str(v))
        j = 0
        for arr in L:
            xy = []
            for corn in arr[0:4]:
                fk5 = corn.toFk5()
                ra, dec = fk5.getRa().asDegrees(), fk5.getDec().asDegrees()
                ras.append(ra)
                decs.append(dec)
                xy.append([ra, dec])
            xy = numpy.array(xy)
            xycen = xy.mean(axis=0)
            poly = Polygon(xy, closed=True, fill=True)
            patches.append(poly)
            colors.append(color)
            ccd_labels.append([xycen[0], xycen[1], arr[4], color])
            if j == 0:
                r = Rectangle((0.0, 0.0), 1, 1, color=cm.jet(color), alpha=0.4)
                patch_proxies.append(r)
                proxy_colors.append(color)
            j += 1
        i += 1

    
    ###########################
    fig = figure.Figure(figsize=(10,10))
    canvas = FigCanvas(fig)
    ax = fig.add_subplot(111)

    ras, decs = [], []
    for k,v in boresights.items():
        ra, dec = v
        c = afwCoord.Fk5Coord(ra, dec)
        cra, cdec = c.getRa().asDegrees(), c.getDec().asDegrees()
        ras.append(cra)
        decs.append(cdec)
        
    if pcen:
        xlo, xhi = pcen.getRa().asDegrees() - rad, pcen.getRa().asDegrees() + rad
        ylo, yhi = pcen.getDec().asDegrees() - rad, pcen.getDec().asDegrees() + rad
    else:
        xlo, xhi = min(ras),  max(ras)
        ylo, yhi = min(decs), max(decs)
        
    if boresightOnly:
        for k,v in boresights.items():
            ra, dec = v
            c = afwCoord.Fk5Coord(ra, dec)
            cra, cdec = c.getRa().asDegrees(), c.getDec().asDegrees()
            offset = -0.05*numpy.abs(xhi - xlo) #ccdMeanWidth
            ax.text(cra+offset, cdec, k, size='xx-small', horizontalalignment='left', verticalalignment='center')
        ax.scatter(numpy.array(ras), numpy.array(decs))
    else:
        pc = PatchCollection(patches, cmap=cm.jet, alpha=0.4)
        pc.set_array(numpy.array(colors))
        ax.add_collection(pc)

    if showLabels and not boresightOnly:
        for ccd_label in ccd_labels:
            xlab, ylab, lab, clr = ccd_label
            ax.text(xlab, ylab, lab, horizontalalignment='center', verticalalignment='center', color=cm.jet(clr))

    if not boresightOnly:
        if n > 1:
            ax.legend(patch_proxies, patch_labels, prop={'size': "large"})
        else:
            ax.set_title(visits[0])

    if pcen:
        xlo, xhi = pcen.getRa().asDegrees() - rad, pcen.getRa().asDegrees() + rad
        ylo, yhi = pcen.getDec().asDegrees() - rad, pcen.getDec().asDegrees() + rad
    else:
        xlo, xhi = min(ras),  max(ras)
        ylo, yhi = min(decs), max(decs)
        
    buff = 0.1
    dx = buff*(xhi - xlo) + edge
    dy = buff*(yhi - ylo) + edge
    
    ax.set_xlim([xhi + dx, xlo - dx])
    ax.set_ylim([ylo - dy, yhi + dy])

    ax.set_xlabel("R.A. [deg]", fontsize='xx-large')
    ax.set_ylabel("Decl. [deg]", fontsize='xx-large')

    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_size('x-large')
    
    if pcen:
        ax.set_title("Boresights")
    else:
        #ax.set_title("chips")
        pass
    
    fig.savefig(out)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action='version', version="0.0")
    parser.add_argument("visit", type=str, help="Visit(s)")
    parser.add_argument("rerun", type=str, help="Rerun to look in.")
    parser.add_argument("-b", "--boresightOnly", action='store_true', default=False,
                        help="Show only the boresight location")
    parser.add_argument("-e", "--edge", type=float, default=0.9,
                        help="edge width to image border (from boresight)")
    parser.add_argument("-l", "--labels", action='store_true', default=False, help="Show CCD labels")
    parser.add_argument("-m", "--minchips", type=int, default=1, help="Only show if n chips within r")
    parser.add_argument("-p", "--point", type=str, default=None, help="show only within r of this point")
    parser.add_argument("-c", "--ccd", type=str, default="0..103",  help="")
    parser.add_argument("-r", "--root", type=str, default=None,  help="Root directory of the data repo")
    parser.add_argument("-o", "--out",  type=str, default="overlap.png",  help="")
    args = parser.parse_args()

    main(args.visit, args.rerun, args.out, ccd=args.ccd,
         root=args.root, point=args.point, minchips=args.minchips, showLabels=args.labels,
         boresightOnly=args.boresightOnly, edge=args.edge)

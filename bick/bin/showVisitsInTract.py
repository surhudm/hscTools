#!/usr/bin/env python

import sys, os, re, math
import argparse
import numpy
import matplotlib.pyplot as pyplot

import lsst.daf.persistence  as dafPersist
import lsst.afw.cameraGeom   as camGeom
import lsst.afw.coord        as afwCoord
import lsst.afw.geom         as afwGeom
import lsst.afw.image        as afwImage

import hsc.tools.bick.utils       as hscUtil


def bboxToRaDec(bbox, wcs):
    """Get the corners of a BBox and convert them to lists of RA and Dec."""
    corners = []
    for corner in bbox.getCorners():
        p = afwGeom.Point2D(corner.getX(), corner.getY())
        coord = wcs.pixelToSky(p).toIcrs()
        corners.append([coord.getRa().asDegrees(), coord.getDec().asDegrees()])
    ra, dec = zip(*corners)
    return ra, dec

def percent(values, p=0.5):
    """Return a value a faction of the way between the min and max values in a list."""
    m = min(values)
    interval = max(values) - m
    return m + p*interval




def main(rootDir, tracts, visitsIn=None, ccds=None, patches=None,
         showPatch=False, showTract=False, showFootprint=False, filter=None, out=None):

    butler = dafPersist.Butler(rootDir)

    if not patches:
        patches = ["%d,%d"%(ix,iy) for ix in range(11) for iy in range(11)]


    #####################################
    # get the visits and ccds
    #####################################
        
    visits = set()
    ccdsInVisits = dict()
    files = {}
    for tract in tracts:
        for patch in patches:
            print "getting patch", tract, patch

            try:
                coadd_file = butler.get("deepCoadd_filename", {'tract':tract, "filter":filter, "patch":patch}, immediate=True)[0]
            except:
                continue
            
            
            if re.search("_parent", coadd_file):
                continue
            
            #print coadd_file
            files[patch] = coadd_file
            try:
                coadd = butler.get("deepCoadd", {'tract':tract, "filter":filter, "patch":patch})
                coaddIn = coadd.getInfo().getCoaddInputs()
            except:
                continue
            ccdInputs = coaddIn.ccds

            for v, ccd in zip(ccdInputs.get("visit"), ccdInputs.get("ccd")):
                if v not in ccdsInVisits:
                    ccdsInVisits[v] = []
                ccdsInVisits[v].append(ccd)

            for v in coaddIn.visits:
                visits.add(int(v.getId()))

    if visitsIn:
        visits = visits & set(visitsIn)
        
    nv = len(visits)

    nvals = {1: (1,1), 2: (2,1), 3: (3,1), 4: (2,2)}
    sizes = {1: (7,7), 2: (10, 5), 3: (12,4), 4: (8,8)}
    if nv in nvals:
        nx, ny = nvals[nv]
        figsize = sizes[nv]
    else:
        nx = int(math.floor(numpy.sqrt(nv)))+1
        ny = nv/nx
        if ny == 0:
            ny = 1
        if nv%nx > 0:
            ny += 1
        if nv < 20:
            figsize = (10, 8)
        else:
            figsize = (2*nx, 2*ny)

    print ny, nx
    fig, axes = pyplot.subplots(ny, nx, squeeze=False, sharex=True, sharey=True, figsize=figsize)
    
    ########################
    ###  draw the CCDs
    ########################
    ras, decs = [], []
    for i_v, visit in enumerate(sorted(visits)):
        print i_v, visit, 
        allCcds = [camGeom.cast_Ccd(ccd) for ccd in camGeom.cast_Raft(butler.get("camera")[0])]
        ccdSet = ccds if ccds else set(ccdsInVisits[visit])
        ccdList = [c for c in allCcds if c.getId().getSerial() in ccdSet]

        color = ('r', 'b', 'c', 'g', 'm') #[0] #[i_v%5]
        i_y, i_x = i_v/nx, i_v%nx
        
        ras = []
        decs = []
        for ccd in ccdList:
            bbox = ccd.getAllPixels()
            ccdId = ccd.getId().getSerial()
            nq = ccd.getOrientation().getNQuarter()

            print "  ", ccdId
            # if it's odd (chips 100..103) color it like it's even counterparts
            nq = nq - nq%2
            clr = color[0]

            # show the halves of the camera, and make chip 0 brighter to break the 180deg symmetry
            alpha = 0.4 if (nq == 0 or ccdId == 0 or ccdId == 50) else 0.2

            if ccdId < 104:
                #dataId = {'tract': tracts, 'visit': visit, 'ccd': ccdId}
                dataId = {'visit': visit, 'ccd': ccdId}
                try:
                    wcs = afwImage.makeWcs(butler.get("calexp_md", dataId))
                except:
                    wcs = None

                if wcs:
                    ra, dec = bboxToRaDec(bbox, wcs)
                    ras += ra
                    decs += dec
                    if not showFootprint:
                        axes[i_y][i_x].fill(ra, dec, fill=True, alpha=alpha, color=clr, edgecolor=clr)

            
                        
        axes[i_y][i_x].set_title(str(visit), size='x-small')
            
        if showFootprint:
            perim = convex_hull(zip(ras, decs))
            pra, pdec = zip(*perim)
            axes[i_y][i_x].fill(pra, pdec, fill=False, edgecolor=color[0])
            
        print ""

    
    ##############################
    # draw the skymap
    ###############################
    minRa = min(ras)
    maxRa = max(ras)
    minDec = min(decs)
    maxDec = max(decs)
    tras = []
    tdecs = []
    if showTract:
        skymap = butler.get('deepCoadd_skyMap', {'tract':0})
        tractObjs = [t for t in skymap if (t.getId() in tracts) ]
        
        for tractObj in tractObjs:
            
            print tractObj.getId()
            ra, dec = bboxToRaDec(tractObj.getBBox(), tractObj.getWcs())
            tras += ra
            tdecs += dec

            # add tract,patch grids to all panels
            for i_x in range(nx):
                for i_y in range(ny):
                    axes[i_y][i_x].fill(ra, dec, fill=False, edgecolor='k', lw=1, linestyle='dashed')
                    axes[i_y][i_x].text(percent(ra, 0.03), percent(dec, 0.95), str(tractObj.getId()),
                                        fontsize=6, horizontalalignment='right', verticalalignment='top')
                    if showPatch:
                        for patch in tractObj:
                            ra, dec = bboxToRaDec(patch.getInnerBBox(), tractObj.getWcs())
                            axes[i_y][i_x].fill(ra, dec, fill=False, edgecolor='k', lw=1, linestyle='dashed')
                            if min(ra) > maxRa or max(ra) < minRa or min(dec) > maxDec or max(dec) < minDec:
                                continue
                            patch_str = "%d,%d" % patch.getIndex()
                            axes[i_y][i_x].text(percent(ra), percent(dec, 0.9), patch_str, 
                                                fontsize=6, horizontalalignment='center',
                                                verticalalignment='top')



    ######################
    # final plot stuff
    ######################

    buff = 0.2
    if showTract:
        minRa = max(min(tras),min(ras))
        maxRa = min(max(tras),max(ras))
        minDec = max(min(tdecs), min(decs))
        maxDec = min(max(tdecs), max(decs))
        xlim = maxRa+buff, minRa-buff
        ylim = minDec-buff, maxDec+buff
    else:
        xlim = max(ras)+buff, min(ras)-buff
        ylim = min(decs)-buff, max(decs)+buff


        
    for i_x in range(nx):
        for i_y in range(ny):
            axes[i_y][i_x].set_xlim(xlim)
            axes[i_y][i_x].set_ylim(ylim)
            for tic in axes[i_y][i_x].get_xticklabels() + axes[i_y][i_x].get_yticklabels():
                tic.set_size("x-small")
            for tic in axes[i_y][i_x].get_xticklabels():
                tic.set_rotation(33.0)
                        
    if rootDir[-1] == "/":
        rootDir = rootDir[:-1]
    rerun = os.path.basename(rootDir)
    tractsort = sorted(tracts)
    if len(tracts) == 1:
        tract_str = str(tracts[0])
    else:
        tract_str = str(tractsort[0]) + "-" + str(tractsort[-1])
    fig.suptitle("Tract %s inputs"%(tract_str))
    outfile = out or "coaddIn-%s-%s.png"%(rerun, tract_str)
    fig.savefig(outfile)
    

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root", help="Root directory of data repository")
    parser.add_argument("tract", help="Tract to show")
    parser.add_argument("-v", "--visits", help="specify Visits")
    parser.add_argument("-c", "--ccds", help="specify CCDs")
    parser.add_argument("-p", "--patches", help="specify patches")
    parser.add_argument("-T", "--showTract", action='store_true', default=False,
                        help="Show the tract boundaries")
    parser.add_argument("-P", "--showPatch", action='store_true', default=False,
                        help="Show the patch boundaries")
    parser.add_argument("-F", '--showFootprint', action='store_true', default=False,
                        help="Show the footprint of the visits instead of the CCDs")
    parser.add_argument("-o", "--out", default=None, help="Output filename")
    parser.add_argument("--filter", default='HSC-I', help="Specify filter for coadd")
    args = parser.parse_args()
        
    main(args.root,
         hscUtil.idSplit(args.tract),
         ccds=hscUtil.idSplit(args.ccds),
         visitsIn=hscUtil.idSplit(args.visits),
         patches=hscUtil.idSplit(args.patches),
         showPatch=args.showPatch, showTract=args.showTract, showFootprint=args.showFootprint,
         filter=args.filter,
         out=args.out)


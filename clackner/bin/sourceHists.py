#!/usr/bin/env python
"""
print the source count histogram for a set of dataIds to stdout.
This includes all sources non-deblended sources and doesn't double count things
in patch/tract overlaps.
The magnitude to use is the psf.mag by default, but another 
magnitude (gauss, sinc, kron) can be specified on the command line.
"""

import sys, os, re
import argparse
import numpy as np

import lsst.daf.persistence as dafPer
import lsst.afw.table as afwTable
import lssttools.utils as lsstUtil

def main(rerun, dataIds, root='/lustre/Subaru/SSP', fakes=None,
         magName='psf'):
    
    doCoadd = 'tract' in dataIds[0].keys()
    butler = dafPer.Butler(os.path.join(root, "rerun", rerun))
    mag_bins = np.arange(12,30,0.05)
    mags = np.zeros(len(mag_bins)-1)
    for dataId in dataIds:
        #print dataId
        try:
            sources = butler.get('deepCoadd_src' if doCoadd else 'src',
                                 dataId, immediate=True, 
                                 flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
            cal_md = butler.get('deepCoadd_md' if doCoadd else 'calexp_md', 
                                dataId, immediate=True)
        except:
            #print "  skipping"
            continue

        srcflux = sources.get('flux.'+magName)
        #mask = sources.get('deblend.nchild') == 0
        mask = sources.get('parent')==0
        if doCoadd:
            mask = mask & (sources.get('detect.is-patch-inner')) & \
                   (sources.get('detect.is-tract-inner'))
        srcmag = -2.5*np.log10(srcflux[mask]/cal_md.get('FLUXMAG0'))
        mags += np.histogram(srcmag, bins=mag_bins)[0]    

        
    #print the magnitude bins and the number counts
    for im, m in enumerate(mags):
        print mag_bins[im]+0.025, m


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('rerun')
    parser.add_argument("visits", help="visits or tracts")
    parser.add_argument("ccds", help="CCDS or patches (for coadds)")
    parser.add_argument("-f", "--filt", 
                        default=None, help="filter, only set for tract/patches")
    parser.add_argument("-R", '--root', default="/lustre/Subaru/SSP")

    args = parser.parse_args()
    visits = lsstUtil.idSplit(args.visits)
    ccds = lsstUtil.idSplit(args.ccds)
    if args.filt is None:
         dataIds = [{'visit':v, 'ccd':c} for c in ccds for v in visits]
    else:
        dataIds = [{'tract':t, 'patch':p, 'filter':args.filt} 
                   for p in ccds for t in visits]
    
    main(args.rerun, dataIds, root=args.root, fakes=args.fakes)
    

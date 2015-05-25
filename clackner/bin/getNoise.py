#!/usr/bin/env python
"""
calculates a prints the expected 5sigma limit for PSF point sources.
This assumes the noise in the pixels isn't correlated, which isn't true on the coadds
"""
import sys, os, re
import argparse
import numpy as np

import lsst.daf.persistence as dafPer
import lsst.afw.table as afwTable
import lssttools.utils as lsstUtil


def main(rerun, dataIds, root='/lustre/Subaru/SSP'):
    doCoadd = 'tract' in dataIds[0].keys()
    butler = dafPer.Butler(os.path.join(root, "rerun", rerun))
    
    for dataId in dataIds:
        try:
            calexp = butler.get('deepCoadd' if doCoadd else 'calexp',
                                dataId, immediate=True)
            md = butler.get('deepCoadd_md' if doCoadd else 'calexp_md',
                            dataId, immediate=True)
        except:
            #print "skipping", dataId
            continue
        
        var = calexp.getMaskedImage().getVariance()
        med_var = np.median(var.getArray().ravel())
        

        # the noise effective area of the PSF
        psf = calexp.getPsf().computeImage()
        Sigmaflux2 = np.sum(psf.getArray()**2)
        
        zp = md.get('FLUXMAG0')

        sn_pp = np.sqrt(med_var/Sigmaflux2)
        print dataId['visit' if not doCoadd else 'tract'], 
        print dataId['ccd' if not doCoadd else 'patch'], 
        #prints the median variance, the psf noise effective area,
        #the 5sigma flux limit and the 5sigma flux limit in magnitudes
        print med_var, Sigmaflux2, 5*sn_pp, -2.5*np.log10(5*sn_pp/zp)
        


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
    
    main(args.rerun, dataIds, root=args.root)
    

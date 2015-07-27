#!/usr/bin/env python
"""
Pulls sources from a set of input data (tracts/patches or visits/ccds)
and dumps it a subset of the source table columns into a single fits file

The output columns include the number of children and the calibrated Kron magnitude
as well as ra/dec/id/parent_id.

This is meant as more of a template than a final product.
"""
import sys, os, re
import argparse
import numpy as np

import lsst.daf.persistence as dafPer
import lsst.afw.table as afwTable
import hsc.tools.bick.utils as hscUtil

def main(rerun, dataIds, root='/lustre/Subaru/SSP', fakes=None):
    
    doCoadd = 'tract' in dataIds[0].keys()
    butler = dafPer.Butler(os.path.join(root, "rerun", rerun))
    schema = afwTable.SourceTable.makeMinimalSchema()
    #add additional columns to output table
    schema.addField('kronMag', type=float, doc='kron magnitude')
    schema.addField('nchild', type=int, doc='nchild')
    #schema.addField('parent', type=long, doc='parent id')
    srcList = afwTable.SourceCatalog(schema)

    #iterate over ccds or patches
    for dataId in dataIds:
        print dataId
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
        #mask = sources.get('parent')==0
        mask = np.ones(len(sources), dtype=bool)
        
        #remove duplicates in patch/tract overlap regions
        if doCoadd:
            mask = mask & (sources.get('detect.is-patch-inner')) & \
                   (sources.get('detect.is-tract-inner'))
        srcmag = -2.5*np.log10(srcflux/cal_md.get('FLUXMAG0'))


        #add sources to output table
        for i_s, source in enumerate(sources[mask]):
            s = srcList.addNew()
            s.setId(source.getId())
            s.setRa(source.getRa())
            s.setDec(source.getDec())
            s.set('kronMag', srcmag[mask][i_s])
            s.set('nchild', source.get('deblend.nchild'))
            s.set('parent', source.get('parent'))

    srcList.writeFits('total_cat.fits')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('rerun')
    parser.add_argument("visits", help="visits or tracts")
    parser.add_argument("ccds", help="CCDS or patches (for coadds)")
    parser.add_argument("-f", "--filt", 
                        default=None, help="filter, only set for tract/patches")
    parser.add_argument("-R", '--root', default="/lustre/Subaru/SSP")
    parser.add_argument('-k', '--fakes', default=None, help='fake catalog to exclude')

    args = parser.parse_args()
    visits = hscUtil.idSplit(args.visits)
    ccds = hscUtil.idSplit(args.ccds)
    if args.filt is None:
         dataIds = [{'visit':v, 'ccd':c} for c in ccds for v in visits]
    else:
        dataIds = [{'tract':t, 'patch':p, 'filter':args.filt} 
                   for p in ccds for t in visits]
    
    main(args.rerun, dataIds, root=args.root, fakes=args.fakes)
    

#!/usr/bin/env python
"""
given an input table with ra/dec, returns a culled list of ra/dec that are not near "DETECTED" pixels in an image.
This returns positions that are in blank regions of the image. The pixel radius necessary to consider a region blank is a command-line argument. 
This can operate over a list of visits/ccds or tracts/patchs
The code expects a fits file and will output a new fits file
"""

import sys, os, re
import argparse
import numpy as np
import pyfits

import lsst.daf.persistence as dafPer
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

import hsc.tools.bick.utils as hscUtil

def loadRaDec(data):
    """
    loads ra and dec from a fits file and puts them in a basic schema.
    this grabs only the ra/dec from a fits file and puts them in an lsst.afw.table
    schema that the pipeline can match against
    """
   
    ras = data['ra']
    decs = data['dec']
    try:
        ids = data['ID']
    except:
        ids = range(len(data))
    #turns ra/dec into basic schema and return the schema
    schema = afwTable.SourceTable.makeMinimalSchema()
    table  = afwTable.SourceTable.make(schema)
    scat   = afwTable.SourceCatalog(table)
    for i,(ra,dec,ident) in enumerate(zip(ras,decs,ids)):
        s = scat.addNew()
        s.setId(int(ident))
        s.setRa(float(ra)*afwGeom.degrees)
        s.setDec(float(dec)*afwGeom.degrees)
    return scat
    


def main(rerun, dataIds, fakes, root='/lustre/Subaru/SSP', rad=10):
    
    doCoadd = 'tract' in dataIds[0].keys()
    butler = dafPer.Butler(os.path.join(root, "rerun", rerun))

    #read in fits file, replace with txt file or anything else
    fits = pyfits.open(fakes)
    data = fits[1].data
    radecCat = loadRaDec(data)
    ndata = len(data)
    datamask = np.ones(ndata, dtype=bool)
    ids = data["ID"] if "ID" in data.names else range(len(data))
    idDict = dict(zip(ids, xrange(ndata)))

    for dataId in dataIds:
        print dataId
        try:
            sources = butler.get('deepCoadd_src' if doCoadd else 'src',
                                 dataId, immediate=True, 
                                 flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
            cal_md = butler.get('deepCoadd_md' if doCoadd else 'calexp_md', 
                                dataId, immediate=True)
            calexp = butler.get('deepCoadd' if doCoadd else 'calexp',
                                dataId, immediate=True)
        except:
            print "skipping", dataId
            continue

        if False:
            matches = afwTable.matchRaDec(sources, radecCat, 
                                          3.3*afwGeom.arcseconds)
            for (src, fake, d) in matches:
                datamask[idDict[fake.getId()]] = False
        
        msk = calexp.getMaskedImage().getMask()
        detected = msk.clone()
        detected &= msk.getPlaneBitMask("DETECTED")
        wcs = calexp.getWcs()
        count, good_count = 0, 0
        for i_d, datum in enumerate(radecCat):
            pixCoord = afwGeom.Point2I(wcs.skyToPixel(datum.getCoord()))
            pixBox = afwGeom.BoxI(pixCoord,afwGeom.Extent2I(1,1))
            pixBox.grow(rad)
            pixBox.clip(calexp.getBBox(afwImage.PARENT))
            if pixBox.isEmpty():
                continue
            else:
                count += 1
                subMask = afwImage.MaskU(detected, pixBox, afwImage.PARENT)
                if sum(subMask.getArray().ravel()) != 0:
                    datamask[i_d] = False
                else:
                    good_count += 1
        print count, good_count

    newdata = data[datamask]
    print ndata, len(newdata)
    hdu = pyfits.BinTableHDU(newdata)
    hdu.writeto('blank_sources.fits', clobber=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('rerun')
    parser.add_argument('fakes')
    parser.add_argument("visits", help="visits or tracts")
    parser.add_argument("ccds", help="CCDS or patches (for coadds)")
    parser.add_argument("-f", "--filt", 
                        default=None, help="filter, only set for tract/patches")
    parser.add_argument("-R", '--root', default="/lustre/Subaru/SSP")
    parser.add_argument('-r', '--radius', 
                        type=int, default=20, help='pixel radius to avoid')

    args = parser.parse_args()
    visits = hscUtil.idSplit(args.visits)
    ccds = hscUtil.idSplit(args.ccds)
    if args.filt is None:
         dataIds = [{'visit':v, 'ccd':c} for c in ccds for v in visits]
    else:
        dataIds = [{'tract':t, 'patch':p, 'filter':args.filt} 
                   for p in ccds for t in visits]
    
    main(args.rerun, dataIds, args.fakes, root=args.root, rad=args.radius)
    

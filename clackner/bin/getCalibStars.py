#!/usr/bin/env python
"""
Returns the calib stars for given dataIds
"""
import os, re
import numpy as np

import lsst.afw.table as afwTable
import lsst.daf.persistence as dafPer
import lsst.pex.exceptions

import lssttools.utils as lsstUtil

import astropy.table


def getMag(flux, fluxerr, zeropoint):
    """
    return the magnitude and error
    """
    mag, magerr = -2.5 * np.log10(flux), 2.5/np.log(10.0)*fluxerr/flux
    return (mag.T + zeropoint).T, magerr


def getEllipse(quad):
    """
    returns the semi-major axis, axes ratio and PA for a given quadrupole moment
    """
    e = lsst.afw.geom.ellipses.Axes(quad)
    return e.getA(), e.getB()/e.getA(), e.getTheta() * 180.0/np.pi


def getAstroTable(src, mags=True, nameList=None, zeropoint=0.0):
    """
    returns an astropy table with all the src entries
    if the entries are complex objects, it breaks them down:
      ellipse entries are broken into 
           ellipse_a = semi-major axis
           ellipse_q = axis ratio (always < 1)
           ellipse_theta = rotation of semi-major axis from chip x-axis in degrees 
    if mags is True, returns the magnitudes for all the flux columns
    """
    
    tab = astropy.table.Table()
    allNames = nameList if nameList is not None else src.schema.getNames()
    for name in allNames:
        #for reasons I don't understand a lookup by name is much slower than a lookup by key
        nameKey = src.schema.find(name).getKey()
        try: 
            tab.add_column(astropy.table.Column(name=name,
                                                data=src.get(nameKey)))
        except lsst.pex.exceptions.LsstException:
            if type(src[0].get(nameKey)) is lsst.afw.geom.ellipses.ellipsesLib.Quadrupole:
                reff, q, theta = zip(*[getEllipse(s.get(nameKey)) for s in src])
                tab.add_column(astropy.table.Column(name=name+'_a', data=reff))
                tab.add_column(astropy.table.Column(name=name+'_q', data=q))
                tab.add_column(astropy.table.Column(name=name+'_theta', data=theta))
            elif type(src[0].get(nameKey)) is lsst.afw.coord.coordLib.IcrsCoord:
                x, y= zip(*[(s.get(nameKey).getRa().asDegrees(), 
                             s.get(nameKey).getDec().asDegrees()) for s in src])
                tab.add_column(astropy.table.Column(name=name+'_ra', data=x))
                tab.add_column(astropy.table.Column(name=name+'_dec', data=y))
            else:
                tab.add_column(astropy.table.Column(name=name, 
                                                    data=np.array([s.get(nameKey) for s in src])))
            #report angles in degrees
        if isinstance(src[0].get(nameKey), lsst.afw.geom.Angle):
            tab.remove_column(name)
            tab.add_column(astropy.table.Column(data=[s.get(nameKey).asDegrees()
                                                      for s in src],
                                                dtype=float, name=name))

    if mags:
        #this is a horrible hack, but I don't think we can use the slots, since 
        #not all the fluxes end up in the slots
        for col in tab.colnames:
            if (re.match('^flux\.[a-z]+$', col) or 
                re.match('^flux\.[a-z]+.apcorr$', col) or
                re.match('^cmodel.+flux$', col) or 
                re.match('^cmodel.+flux.apcorr$', col)):
                try:
                    zp = tab['zeropoint']
                except:
                    zp = zeropoint
                zp = 0.0 if re.search('apcorr', col) else zp
                mag, magerr = getMag(tab[col], tab[col+'.err'], zp)
            
                tab.add_column(astropy.table.Column(name=re.sub('flux', 'mag', col),
                                                    data=mag))
                tab.add_column(astropy.table.Column(name=re.sub('flux', 'mag', col+'.err'),
                                                    data=magerr))
                
    return tab



def main(rerun, dataIds, root):
    doCoadd = 'tract' in dataIds[0].keys()
    butler = dafPer.Butler(os.path.join(root, "rerun", rerun))
    
    out_tab = None

    for dataId in dataIds:
        print dataId
        try:
            src = butler.get("deepCoadd_src" if doCoadd else "src", 
                                 dataId, immediate=True)
        except:
            print "no source catalog for: ",repr(dataId)
            continue
        
        isPsfStar = src.get('calib.psf.candidate')
        isPsfReserved = src.get('calib.psf.reserved')
        print  butler.get("calexp", dataId).getCalib().getFluxMag0()

        cutsrc = src[isPsfStar | isPsfReserved]
    
        cutsrctab = getAstroTable(cutsrc, nameList=['id', 'coord', 
                                                    'parent', 'deblend.nchild',
                                                    'flags.pixel.edge',
                                                    'flags.pixel.interpolated.any',
                                                    'shape.sdss',
                                                    'classification.extendedness',
                                                    'calib.psf.used',
                                                    'calib.psf.candidate',
                                                    'calib.psf.reserved',
                                                    'flux.psf', 'flux.gaussian',
                                                    'flux.psf.err', 'flux.gaussian.err'],
                                  zeropoint = 2.5*np.log10(butler.get("calexp", dataId).getCalib().getFluxMag0()[0]))
        if out_tab is None:
            out_tab = astropy.table.Table(cutsrctab, copy=True)
        else:
            out_tab = astropy.table.vstack([out_tab, cutsrctab], join_type='exact')

    out_tab.write("calib_psf.fits", overwrite=True)
        

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("rerun")
    parser.add_argument("visits", help="visits or tracts")
    parser.add_argument("ccds", help="CCDs or patches, separated by ^")
    parser.add_argument("-f", "--filt", 
                        default=None, help="filter, only set for tract/patches")
    parser.add_argument("-R", "--root", default="/lustre/Subaru/SSP")

    args = parser.parse_args()
    visits = lsstUtil.idSplit(args.visits)
    ccds = lsstUtil.idSplit(args.ccds)
    if args.filt is None:
        dataIds = [{'visit':v, 'ccd':c} for v in visits for c in ccds]
    else:
        dataIds = [{'tract':t, 'patch':p, 'filter':args.filt} 
                   for p in ccds for t in visits]
    
    main(args.rerun, dataIds, root=args.root)
    

#!/usr/bin/env python

import hsc.tools.bick.safeTools as safeTools
import lsst.afw.image      as afwImage
import lsst.afw.cameraGeom as camGeom

safeTools.butlerTarget = "src"

class FindTractTask(safeTools.SafeTask):
    """A Task to find the tract associated with a visit/ccd.
    """
    def run(self, dataRef):

        calexp_md = dataRef.get("calexp_md")
        wcs    = afwImage.makeWcs(calexp_md)
        skymap = dataRef.get("deepCoadd_skyMap", immediate=True)

        # all this, just to get the center pixel coordinate
        camera = dataRef.get("camera")
        raft   = camGeom.cast_Raft(camera.findDetector(camGeom.Id(0)))
        detId  = camGeom.Id(calexp_md.get("DET-ID"))
        ccd    = camGeom.cast_Ccd(raft.findDetector(detId))
        size   = ccd.getSize().getPixels(ccd.getPixelSize())
        
        coord  = wcs.pixelToSky(size.getX()/2, size.getY()/2)
        tract  = skymap.findTract(coord).getId()

        d = dataRef.dataId
        print "%-6d %3d  %5d" % (d['visit'], d['ccd'], tract)
        return tract

if __name__ == '__main__':
    FindTractTask.parseAndRun()


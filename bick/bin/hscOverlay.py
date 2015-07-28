#!/usr/bin/env python
# Original filename: bin/hscOverlay.py
#
# Author: Steven Bickerton
# Email: 
# Date: Tue 2014-10-07 11:50:49
# 
# Summary: 
# 

import sys
import os
import re
import math
import argparse
import lsst.daf.persistence          as dafPersist
import hsc.pipe.base.butler          as hscButler
import numpy
import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
import PIL

import hsc.tools.bick.utils as hscUtil

import lsst.afw.cameraGeom          as camGeom
import lsst.afw.cameraGeom.utils    as camGeomUtils



#############################################################
#
# Main body of code
#
#############################################################

def main(figname, camera='hsc', output=None):


    # build a camera
    if camera == 'hsc':
        simdir         = os.environ['OBS_SUBARU_DIR']
        cameraGeomPaf  = os.path.join(simdir, "hsc", "hsc_geom.paf")
        if not os.path.exists(cameraGeomPaf):
            cameraGeomPaf = os.path.join(simdir, "hsc", "description", "hsc_geom.paf")
            if not os.path.exists(cameraGeomPaf):
                raise Exception("Unable to find cameraGeom Policy file: %s" % (cameraGeomPaf))
        cameraGeomPolicy = camGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = camGeomUtils.makeCamera(cameraGeomPolicy)
        
    if not camera:
        raise ValueError("Camera must be ... uhm ... 'hsc'")

    ###########################
    # make RBG arrays with FpaImage objects
    # The FpaImage objects have a method to extract pixels from one CCD
    # So we'll create a dummy image and put the whole user image in it,
    # then we'll copy each CCD into the final image
    ###########################
    bin = 16
    fpa_img = [
        hscUtil.FpaImage(camera, scale=bin),
        hscUtil.FpaImage(camera, scale=bin),
        hscUtil.FpaImage(camera, scale=bin),
    ]
    fpa_dum = [
        hscUtil.FpaImage(camera, scale=bin),
        hscUtil.FpaImage(camera, scale=bin),
        hscUtil.FpaImage(camera, scale=bin),
    ]
    ny, nx = fpa_img[0].image.shape

    ###########################
    # Load the image the user provided.
    # figure out how much to scale things so it fits in our FpaImage
    ###########################
    img0 = PIL.Image.open(figname)
    h, w = img0.size
    if h >= w:
        scale = 1.0*ny/h
    else:
        scale = 1.0*nx/w
    img = numpy.array(img0.resize((int(scale*h),int(scale*w)), PIL.Image.ANTIALIAS))
    subw = int(scale*w)
    subh = int(scale*h)
    x0 = int(0.5*(nx - subw))
    y0 = int(0.5*(ny - subh))

    fpa_dum[0].image += 100
    for i in 0,1,2:
        fpa_dum[i].image[y0:y0+subh,x0:x0+subw] = img[:,:,i]

    ###########################
    # Now copy each CCD from the dummy to the original
    ###########################
    for i in 0,1,2: #, 1, 2:
        for i_ccd in range(104):
            img2 = fpa_dum[i].getPixels(i_ccd)
            print i_ccd, img2.shape
            fpa_img[i].insert(i_ccd, img2)

    ###########################
    # convert this to an RBG image
    # we don't need to worry about uint8 with range 256, matplotlib will handle normalization.
    ###########################
    color_img = numpy.ones((ny, nx, 3), dtype=numpy.uint8)
    for i in range(3):
        color_img[:,:,i] = fpa_img[i].image

        
    ###########################
    # plot it
    ###########################
    fig = figure.Figure()
    canvas = FigCanvas(fig)
    ax = fig.add_subplot(111)
    ax.imshow(color_img)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    fig.savefig(output)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action='version', version="0.0")
    parser.add_argument("figure", type=str, help="")
    parser.add_argument("-c", "--camera", default='hsc', choices=("hsc"), help="Which camera?")
    parser.add_argument("-o", "--output", default=None, help="")
    args = parser.parse_args()

    if not args.output:
        fname, ext = os.path.splitext(args.figure)
        args.output = fname + "-hsc" + ext
    
    main(args.figure, output=args.output, camera=args.camera)

#!/usr/bin/env python

import sys
import optparse
import lsst.afw.image        as afwImage
import lsst.afw.display.ds9  as ds9

"""
%prog [options] arg
"""


def main(infiles, frame=1, title="", scale="zscale", zoom="to fit", trans=60):

    settings = {'scale':scale, 'zoom': zoom, 'mask' : 'transparency %d' %(trans)}

    for infile in infiles:
        t = title
        if not title:
            t = infile
        mimg = afwImage.MaskedImageF(infile)
        ds9.mtv(mimg, frame=frame, title=t, settings=settings)
        frame += 1
    
if __name__ == '__main__':

    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("-f", "--frame", dest="frame", type=int,
                      default=1, help="Frame (default=%default)")
    parser.add_option("-s", "--scale",
                      default="zscale", help="Scale (default=%default)")
    parser.add_option("-t", "--title", dest="title",
                      default="", help="Title (default=%default)")
    parser.add_option("-T", "--trans",
                      default=60, help="Transparency (default=%default)")
    parser.add_option("-z", "--zoom",
                      default="to fit", help="Zoom (default=%default)")
    
    opts, args = parser.parse_args()


    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    fitsfiles = args

    main(fitsfiles, frame=opts.frame, title=opts.title, scale=opts.scale, zoom=opts.zoom, trans=opts.trans)

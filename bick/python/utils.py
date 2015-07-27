
import sys, os, re
import numpy

#import matplotlib.figure as figure
#from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
import matplotlib.patches   as patches

import lsst.afw.geom        as afwGeom
import lsst.afw.cameraGeom  as camGeom
import lsst.daf.persistence as dafPersist

def idSplit(id):
    if id is None:
        return id
    ids = []
    for r in id.split("^"):
        m = re.match(r"^(\d+)\.\.(\d+):?(\d+)?$", r)
        if m:
            limits = [int(v) if v else 1 for v in m.groups()]
            limits[1] += 1
            ids += range(*limits)
        else:
            ids.append(r if r.count(",") else int(r))
    return ids


def numSplit(s, discrete_delim=',', range_delim='-'):
    nums = []
    chunks = s.split(discrete_delim)
    for c in chunks:
        if re.search(range_delim, c):
            n = range(*map(int, c.split(range_delim)))
        else:
            n = [int(c)]
        nums += n
    return nums


def rebin(img, b, mean=True):
    ny, nx = img.shape
    tmp = numpy.zeros((ny/b, nx/b))
    mx = nx % b
    my = ny % b
    for i in range(b):
        for j in range(b):
            tmp += img[i:ny-my:b,j:nx-mx:b]
    if mean:
        tmp /= b*b
    return tmp


    
def getButler(rerun, root=None):
    
    ############################
    dirOptions = []
    if root and os.path.join(root, rerun):
        dirOptions.append(os.path.join(root, rerun))
    else:
        if "DATA_PATH" not in os.environ:
            raise RuntimeError("Must set DATA_PATH or specify --root=/path/to/data")
        dirs = os.environ["DATA_PATH"].split(":")
        for d in dirs:
            r = os.path.join(d, rerun)
            if os.path.exists(r):
                dirOptions.append(r)
    if len(dirOptions) == 0:
        raise RuntimeError("Unable to find rerun " + rerun)
        
    run_index = 0
    if len(dirOptions) > 1:
        data_select = os.environ.get("DATA_SELECT", "query")

        if data_select == "query":
            msg = "Rerun " + rerun + " is ambiguous\n" + '\n'.join(["(%d) %s"%(i, d) for i,d in enumerate(dirOptions)])
            print msg
            run_index = int(raw_input("Which one do I run? "))
        if data_select == 'first':
            run_index = 0
        if data_select == 'throw':
            raise RuntimeError(msg)

        
    if root is None:
        root = os.path.dirname(dirOptions[run_index])
        
    if os.path.exists(os.path.join(root, rerun, '_parent')):
        rootDir = dirOptions[run_index]
    elif os.path.exists(os.path.join(root, '_mapper')):
        rootDir = root
    print "# Getting butler using root = ", rootDir
    butler = dafPersist.Butler(rootDir)

    return butler


    
def histeq(inData):

    inShape = inData.shape
    data = inData.flatten()
    n = len(data)
    temp = numpy.argsort(data)
    ranks = numpy.empty(n, int)
    ranks[temp] = numpy.arange(n)
    histeq = ranks.reshape(inShape)
    return histeq

def zscale(img, contrast=0.25, samples=500):
    ravel = img.ravel()
    if len(ravel) > samples:
        imsort = numpy.sort(numpy.random.choice(ravel, size=samples))
    else:
        imsort = numpy.sort(ravel)

    n = len(imsort)
    idx = numpy.arange(n)

    med = imsort[n/2]
    w = 0.25
    i_lo, i_hi = int((0.5-w)*n), int((0.5+w)*n)
    p = numpy.polyfit(idx[i_lo:i_hi], imsort[i_lo:i_hi], 1)
    slope, intercept = p

    z1 = med - (slope/contrast)*(n/2-n*w)
    z2 = med + (slope/contrast)*(n/2-n*w)

    return z1, z2    
    

class FpaImage(object):

    def __init__(self, camera, scale=16.0):
        self.camera  = camera
        self.detectors = {}
        self.amps      = {}

        self.xy0 = {}
        self.dims = {}
        self.used = set()
        self.nQuarter = {}
        
        l_cam, b_cam, r_cam, t_cam = 1.0e6, 1.0e6, -1.0e6, -1.0e6
        for r in self.camera:
            raft = camGeom.cast_Raft(r)
            for c in raft:
                ccd = camGeom.cast_Ccd(c)
                serial = ccd.getId().getSerial()
                self.detectors[serial] = ccd

                cc       = ccd.getCenter().getPixels(ccd.getPixelSize())
                cxc, cyc = cc.getX(), cc.getY()
                cbbox    = ccd.getAllPixels(True)
                cwidth   = cbbox.getMaxX() - cbbox.getMinX()
                cheight  = cbbox.getMaxY() - cbbox.getMinY()

                x0, y0 = int(cxc - cwidth/2.0), int(cyc - cheight/2.0)

                self.xy0[serial] = (x0, y0)
                self.dims[serial] = (cwidth/scale, cheight/scale)
                self.nQuarter[serial] = ccd.getOrientation().getNQuarter()
                
                l_cam = min(l_cam, x0)
                b_cam = min(b_cam, y0)
                r_cam = max(r_cam, cxc + cwidth/2)
                t_cam = max(t_cam, cyc + cheight/2)

        for k,v in self.xy0.items():
            self.xy0[k] = [int(x) for x in ((v[0] - l_cam)/scale, (v[1] - b_cam)/scale)]
            
        self.camera_width  = r_cam - l_cam
        self.camera_height = t_cam - b_cam

        extent = afwGeom.Extent2D(self.camera_width, self.camera_height)
        self.camera_bbox   = afwGeom.Box2D(afwGeom.Point2D(0.0, 0.0), extent)

        self.image = numpy.zeros((self.camera_height/scale, self.camera_width/scale))

    def insert(self, ccd, img, bin=None):
        x0, y0 = self.xy0[ccd]
        w, h = self.dims[ccd]
        nq = self.nQuarter[ccd]

        ih, iw = img.shape
        if not bin:
            bin = int(ih*1.0/h)
        img_bin = rebin(img, bin)[:h,:w]
        self.image[y0:y0+h,x0:x0+w] = img_bin
        self.used.add(ccd)

    def getPixels(self, ccd):
        x0, y0 = self.xy0[ccd]
        w, h = self.dims[ccd]
        nq = self.nQuarter[ccd]
        return self.image[y0:y0+h,x0:x0+w]
        
    def getValidPixels(self):
        tmp = []
        for ccd in self.used:
            x0, y0 = self.xy0[ccd]
            w, h = self.dims[ccd]
            nq = self.nQuarter[ccd]
            tmp.append(self.image[y0:y0+h,x0:x0+w])
        return numpy.array(tmp)

        
class FpaFigure(object):

    def __init__(self, figure, camera, rect=(0.08, 0.08, 0.84, 0.84)):
        self.figure  = figure
        self.camera  = camera
        self.rect    = rect
        self.detectors = {}
        self.amps      = {}
        
        # use same rect as figure.add_axes(rect) to avoid confusion
        l, b, w, h = self.rect
        self.fig_left   = l
        self.fig_bottom = b
        self.fig_width  = w
        self.fig_height = h

        l_cam, b_cam, r_cam, t_cam = 1.0e6, 1.0e6, -1.0e6, -1.0e6
        for r in self.camera:
            raft = camGeom.cast_Raft(r)
            for c in raft:
                ccd = camGeom.cast_Ccd(c)
                serial = ccd.getId().getSerial()
                self.detectors[serial] = ccd

                cc       = ccd.getCenter().getPixels(ccd.getPixelSize())
                cxc, cyc = cc.getX(), cc.getY()
                cbbox    = ccd.getAllPixels(True)
                cwidth   = cbbox.getMaxX() - cbbox.getMinX()
                cheight  = cbbox.getMaxY() - cbbox.getMinY()

                l_cam = min(l_cam, cxc - cwidth/2)
                b_cam = min(b_cam, cyc - cheight/2)
                r_cam = max(r_cam, cxc + cwidth/2)
                t_cam = max(t_cam, cyc + cheight/2)
                
                self.amps[serial] = []
                for a in ccd:
                    amp = camGeom.cast_Amp(a)
                    self.amps[serial].append(amp)


        self.camera_width  = r_cam - l_cam
        self.camera_height = t_cam - b_cam

        extent = afwGeom.Extent2D(self.camera_width, self.camera_height)
        self.camera_bbox   = afwGeom.Box2D(afwGeom.Point2D(0.0, 0.0), extent)

        #print "Camera", self.camera_bbox
        
                    
        self.axes = {}
                
    def getAxes(self, serial):
        ccd = self.detectors[serial]
        raft = camGeom.cast_Raft(ccd.getParent())

        # NOTE: all ccd coords are w.r.t. the *center* of the raft, not its LLC
        if False:
            rc        = raft.getCenter().getPixels(raft.getPixelSize())
            rxc, ryc  = rc.getX(), rc.getY()
            if rxc == 0 or ryc == 0:
                rbbox = raft.getAllPixels(True)
                rxc   = rbbox.getWidth()/2
                ryc   = rbbox.getHeight()/2
        else:
            rxc = self.camera_width/2
            ryc = self.camera_height/2

        #label = ccd.getId().getName()

        cc       = ccd.getCenter().getPixels(ccd.getPixelSize())
        cxc, cyc = cc.getX(), cc.getY()
        cbbox    = ccd.getAllPixels(True)
        cwidth   = cbbox.getMaxX() - cbbox.getMinX()
        cheight  = cbbox.getMaxY() - cbbox.getMinY()

        #print rxc, ryc, cxc, cyc
        
        cx0      = rxc + cxc - cwidth/2
        cy0      = ryc + cyc - cheight/2
        cx1      = cx0 + cwidth
        cy1      = cy0 + cheight

        # rescale according to rect
        sx = self.camera_width/self.fig_width
        sy = self.camera_height/self.fig_height
        l, w = self.fig_left   + cx0/sx, cwidth/sx
        b, h = self.fig_bottom + cy0/sy, cheight/sy
        rect = (l, b, w, h)

        ax = self.figure.add_axes(rect)

        # kill the ticks
        for tick in ax.get_xticklines() + ax.get_yticklines() + ax.get_yticklabels() + ax.get_xticklabels():
            tick.set_visible(False)

        self.axes[serial] = ax
        return ax

    def addLabel(self, serial, labels, color='#00ff00', size=8):

        ax = self.axes.get(serial, self.getAxes(serial))

        # labels
        nQuarter = self.detectors[serial].getOrientation().getNQuarter()
        offset = 0.15
        if nQuarter % 2 == 1:
            offset = 0.3
        
        for i in range(len(labels)):
            ax.text(0.5, 0.9 - i*offset, labels[i], color=color, fontsize=size, 
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    def highlightAmp(self, serial, amp, color='red'):

        ax = self.axes.get(serial, self.getAxes(serial))
        ccd_bb = self.detectors[serial].getAllPixels()
        amp_bb = self.amps[serial][amp].getAllPixels()
        #xmin, xmax = amp_bb.getMinX(), amp_bb.getMaxX()
        #ymin, ymax = amp_bb.getMinY(), amp_bb.getMaxY()
        buff = 0.02
        xmin = 1.0*amp_bb.getMinX()/ccd_bb.getWidth()+buff
        xmax = 1.0*amp_bb.getMaxX()/ccd_bb.getWidth()-buff
        ymin = 1.0*amp_bb.getMinY()/ccd_bb.getHeight()+buff
        ymax = 1.0*amp_bb.getMaxY()/ccd_bb.getHeight()-buff

        rect = patches.Rectangle((xmin,ymin), xmax-xmin, ymax-ymin, edgecolor=color, facecolor='none',
                                 transform=ax.transAxes)
        ax.add_patch(rect)

        
def showCcdAxes(fig, camera, rect, serials):

    fpa_fig = FpaFigure(fig, camera, rect=rect)
    for s in fpa_fig.detectors.keys():
        ax = fpa_fig.getAxes(s)
        if s in serials:
            ax.patch.set_facecolor('k')


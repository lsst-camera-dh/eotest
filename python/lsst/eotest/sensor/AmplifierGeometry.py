"""
Code to encapsulate amplifier geometry as expressed in NOAO image
section keywords DETSEC, DATASEC, DETSIZE.
"""
import astropy.io.fits as fits

amp_loc = {}
amp_loc['E2V'] = dict([(amp, -1) for amp in range(1, 9)] +
                      [(amp, 1) for amp in range(9, 17)])
amp_loc['ITL'] = dict([(amp, -1) for amp in range(1, 9)] +
                      [(amp, -1) for amp in range(9, 17)])

def parse_geom_kwd(value):
    geom = {}
    data = value[1:-1].split(',')
    xmin, xmax = (int(x) for x in data[0].split(':'))
    geom['xmin'] = xmin
    geom['xmax'] = xmax
    ymin, ymax = (int(y) for y in data[1].split(':'))
    geom['ymin'] = ymin
    geom['ymax'] = ymax
    return geom

def makeAmplifierGeometry(infile):
    """
    Make an AmplifierGeometry object from an input FITS file.
    """
    foo = fits.open(infile)
    detsize = parse_geom_kwd(foo[0].header['DETSIZE'])
    datasec = parse_geom_kwd(foo[1].header['DATASEC'])
    prescan = datasec['xmin'] - 1
    nx = datasec['xmax'] - prescan
    ny = datasec['ymax'] - datasec['ymin'] + 1
    # kludge to infer amplifier node locations for two vendors
    detsec = parse_geom_kwd(foo[9].header['DETSEC'])
    if detsec['xmin'] < detsec['xmax']:
        vendor = 'E2V'
    else:
        vendor = 'ITL'
    myAmpGeom = AmplifierGeometry(prescan=prescan, nx=nx, ny=ny,
                                  detxsize=detsize['xmax'],
                                  detysize=detsize['ymax'],
                                  amp_loc=amp_loc[vendor])
    myAmpGeom.compute_geometry(fitsfile=infile)
    return myAmpGeom
    
class AmplifierGeometry(dict):
    nsegx, nsegy = 8, 2
    def __init__(self, prescan=10, nx=512, ny=2002,
                 detxsize=4336, detysize=4044, amp_loc=amp_loc['E2V']):
        super(AmplifierGeometry, self).__init__()
        self.prescan_width = prescan
        self.nx = nx
        self.ny = ny
        self.amp_loc = amp_loc
        self.compute_geometry(detxsize=detxsize, detysize=detysize)
    def compute_geometry(self, **kwds):
        if kwds.has_key('fitsfile'):
            # Compute geometry by inferring DETSIZE from NAXIS[12] in
            # first image extension of specified FITS file.
            foo = fits.open(kwds['fitsfile'])
            self.naxis1 = foo[1].header['NAXIS1']
            self.naxis2 = foo[1].header['NAXIS2']
            detxsize = self.naxis1*self.nsegx
            detysize = self.naxis2*self.nsegy
        else:
            # Compute geometry using supplied detxsize, detysize
            detxsize = kwds['detxsize']
            detysize = kwds['detysize']
            self.naxis1 = detxsize/self.nsegx
            self.naxis2 = detysize/self.nsegy
        self.DETSIZE = '[1:%i,1:%i]' % (detxsize, detysize)
        self.serial_overscan_width = self.naxis1 - self.nx - self.prescan_width
        try:
            self._make_bboxes()
        except ImportError:
            pass
        for amp in range(1, self.nsegx*self.nsegy + 1):
            self[amp] = self._segment_geometry(amp)
    def _make_bboxes(self):
        import lsst.afw.geom as afwGeom
        self.full_segment = \
            afwGeom.Box2I(afwGeom.Point2I(0, 0),
                          afwGeom.Point2I(self.naxis1 - 1, self.naxis2 - 1))
        self.prescan = \
            afwGeom.Box2I(afwGeom.Point2I(0, 0),
                          afwGeom.Point2I(self.prescan_width - 1,
                                          self.naxis2 - 1))
        self.imaging = \
            afwGeom.Box2I(afwGeom.Point2I(self.prescan_width, 0),
                          afwGeom.Point2I(self.nx + self.prescan_width - 1,
                                          self.ny - 1 ))
        self.serial_overscan = \
            afwGeom.Box2I(afwGeom.Point2I(self.nx + self.prescan_width, 0),
                          afwGeom.Point2I(self.naxis1 - 1, self.naxis2 - 1))
        self.parallel_overscan = \
            afwGeom.Box2I(afwGeom.Point2I(self.prescan_width, self.ny),
                          afwGeom.Point2I(self.prescan_width + self.nx,
                                          self.naxis2 - 1))
    def _segment_geometry(self, amp):
        results = dict()
        results['DETSIZE'] = '[1:%i,1:%i]' % (self.nx*self.nsegx, 
                                              self.ny*self.nsegy)
        results['DATASEC'] = \
            '[%i:%i,%i:%i]' % (self.prescan_width + 1,
                               self.naxis1 - self.serial_overscan_width,
                               1, self.ny)
        results['DETSEC'] = self._detsec(amp)
        results['BIASSEC'] = \
            '[%i:%i,%i:%i]' % (self.prescan_width + self.nx + 1,
                               self.naxis1,
                               1, self.ny)
        return results
    def _detsec(self, amp):
        namps = self.nsegx*self.nsegy
        if amp <= self.nsegx:
            x1 = (amp - 1)*self.nx + 1
            x2 = amp*self.nx
            y1, y2 = 1, self.ny
        else: 
            # Amps in "top half" of CCD, where the ordering of amps 9
            # to 16 is right-to-left.
            x1 = (namps - amp)*self.nx + 1
            x2 = (namps - amp + 1)*self.nx
            # Flip in y for top half of sensor.
            y1, y2 = 2*self.ny, self.ny + 1  
        if self.amp_loc[amp] < 0:
            # Flip since the output node is on the right side of segment.
            x1, x2 = x2, x1
        return '[%i:%i,%i:%i]' % (x1, x2, y1, y2)
    def __eq__(self, other):
        for key in self.__dict__.keys():
            if getattr(self, key) != getattr(other, key):
                return False
        return True

if __name__ == '__main__':
    e2v = AmplifierGeometry()
    itl = AmplifierGeometry(prescan=3, nx=509, ny=2000, amp_loc=amp_loc['ITL'])

    print e2v.full_segment
    print e2v.prescan
    print e2v.imaging
    print e2v.serial_overscan
    print e2v.parallel_overscan

    for amp in range(1, 17):
        print amp, e2v[amp]['DETSEC'], itl[amp]['DETSEC']
        print amp, e2v[amp]['DATASEC'], itl[amp]['DATASEC']
        print 
    infile = '/u/gl/jchiang/ki18/LSST/SensorTests/eotest/0.0.0.6/work/sensorData/000-00/fe55/debug/000-00_fe55_bias_00_debug.fits'
    geom = makeAmplifierGeometry(infile)

    print geom == e2v
    print geom == itl

    print geom != e2v
    print geom != itl

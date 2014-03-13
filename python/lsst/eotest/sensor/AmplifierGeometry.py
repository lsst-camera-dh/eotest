import lsst.afw.geom as afwGeom

e2v_amp_loc = dict([(amp, -1) for amp in range(1, 9)] +
                   [(amp, 1) for amp in range(9, 17)])
itl_amp_loc = dict([(amp, -1) for amp in range(1, 9)] +
                   [(amp, -1) for amp in range(9, 17)])

class AmplifierGeometry(dict):
    nsegx, nsegy = 8, 2
    def __init__(self, prescan=10, nx=512, ny=2002,
                 detxsize=4336, detysize=4044, amp_loc=e2v_amp_loc):
        super(AmplifierGeometry, self).__init__()
        self.nx = nx
        self.ny = ny
        self.naxis1 = detxsize/self.nsegx
        self.naxis2 = detysize/self.nsegy
        self.amp_loc = amp_loc
        self.full_segment = \
            afwGeom.Box2I(afwGeom.Point2I(0, 0),
                          afwGeom.Point2I(self.naxis1 - 1, self.naxis2 - 1))
        self.prescan = \
            afwGeom.Box2I(afwGeom.Point2I(0, 0),
                          afwGeom.Point2I(prescan - 1, self.naxis2 - 1))
        self.imaging = \
            afwGeom.Box2I(afwGeom.Point2I(prescan, 0),
                          afwGeom.Point2I(self.nx + prescan - 1, self.ny - 1 ))
        self.serial_overscan = \
            afwGeom.Box2I(afwGeom.Point2I(self.nx + prescan, 0),
                          afwGeom.Point2I(self.naxis1 - 1, self.naxis2 - 1))
        self.parallel_overscan = \
            afwGeom.Box2I(afwGeom.Point2I(prescan, ny),
                          afwGeom.Point2I(prescan + nx, self.naxis2 - 1))
        self.DETSIZE = '[1:%i,1:%i]' % (detxsize, detysize)
        for amp in range(1, self.nsegx*self.nsegy + 1):
            self[amp] = self._segment_geometry(amp)
    def _segment_geometry(self, amp):
        results = dict()
        results['DETSIZE'] = '[1:%i,1:%i]' % (self.nx*self.nsegx, 
                                              self.ny*self.nsegy)
        results['DATASEC'] = \
            '[%i:%i,%i:%i]' % (self.prescan.getWidth() + 1,
                               self.naxis1 - self.serial_overscan.getWidth(),
                               1, self.ny)
        results['DETSEC'] = self._detsec(amp)
        return results
    def _detsec(self, amp):
        namps = self.nsegx*self.nsegy
        if self.amp_loc[amp] < 0:
            x1 = self.nx*amp
            x2 = x1 - self.nx + 1
        else:
            x1 = self.nx*(namps - amp) + 1
            x2 = x1 + self.nx - 1 
        if amp < self.nsegx + 1:
            y1, y2 = 1, self.ny
        else:
            y1, y2 = 2*self.ny, self.ny + 1
        return '[%i:%i,%i:%i]' % (x1, x2, y1, y2)

if __name__ == '__main__':
    e2v = AmplifierGeometry()
    itl = AmplifierGeometry(prescan=3, nx=509, ny=2000, amp_loc=itl_amp_loc)

    print e2v.full_segment
    print e2v.prescan
    print e2v.imaging
    print e2v.serial_overscan
    print e2v.parallel_overscan

    for amp in range(1, 17):
        print amp, e2v[amp]['DETSEC'], itl[amp]['DETSEC']
        print amp, e2v[amp]['DATASEC'], itl[amp]['DATASEC']
        print 


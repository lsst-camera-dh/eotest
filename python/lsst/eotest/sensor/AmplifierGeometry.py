class AmplifierGeometry(dict):
    nsegx, nsegy = 8, 2
    def __init__(self, prescan=10, serial_overscan=20, parallel_overscan=20,
                 detxsize=4336, detysize=4044):
        super(AmplifierGeometry, self).__init__()
        self.prescan = prescan
        self.serial_overscan = serial_overscan
        self.parallel_overscan = parallel_overscan
        self.naxis1 = detxsize/self.nsegx
        self.naxis2 = detysize/self.nsegy
        self.nx = self.naxis1 - self.prescan - self.serial_overscan
        self.ny = self.naxis2 - self.parallel_overscan
        self.DETSIZE = '[1:%i,1:%i]' % (detxsize, detysize)
        for amp in range(1, self.nsegx*self.nsegy+1):
            self[amp] = self._segment_geometry(amp)
    def _segment_geometry(self, amp):
        results = dict()
        results['DETSIZE'] = '[1:%i,1:%i]' % (self.nx*self.nsegx, self.ny*self.nsegy)
        results['DATASEC'] = '[%i:%i,%i:%i]' % (self.prescan + 1,
                                                self.naxis1 - self.serial_overscan,
                                                1, self.naxis2 - self.parallel_overscan)
        results['DETSEC'] = self._detsec(amp)
        return results
    def _detsec(self, amp):
        namps = self.nsegx*self.nsegy
        if amp < self.nsegx + 1:
            x1 = self.nx*amp
            x2 = x1 - self.nx + 1
            y1, y2 = 1, self.ny
        else:
            x1 = self.nx*(namps - amp) + 1
            x2 = x1 + self.nx - 1 
            y1, y2 = 2*self.ny, self.ny + 1
        return '[%i:%i,%i:%i]' % (x1, x2, y1, y2)

if __name__ == '__main__':
    e2v = AmplifierGeometry()
    itl = AmplifierGeometry(prescan=3, detxsize=4400, detysize=4040)
    for amp in range(1, 17):
        print amp, e2v[amp]['DETSEC'], itl[amp]['DETSEC']
        print amp, e2v[amp]['DATASEC'], itl[amp]['DATASEC']
        print 

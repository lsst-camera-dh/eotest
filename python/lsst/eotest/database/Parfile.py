"""
@brief Class to parse named parameters in a text file as key/value
pairs and return the output as a dict subclass, performing conversions
to int and float as appropriate.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os


class Parfile(dict):
    def __init__(self, filename, fixed_keys=True):
        """
        fixed_keys=True indicates that the input file has all of the
        valid keys already specified and that new keys cannot be
        added.
        """
        self.filename = os.path.abspath(filename)
        self.fixed_keys = fixed_keys
        self.keylist = []
        try:
            self._readfile()
        except IOError:
            if fixed_keys:
                raise
            pass

    def _readfile(self):
        for line in open(self.filename):
            if line.find('#') == 0:
                continue
            if line.find(" #") != -1 or line.find("\t#") != -1:
                data = '#'.join(line.split('#')[:-1])
            else:
                data = line
            key, value = [x.strip() for x in data.split("=")]
            self._addkey(key)
            try:
                try:
                    self[key.strip()] = int(value.strip())
                except:
                    self[key.strip()] = float(value.strip())
            except ValueError:
                self[key.strip()] = value.strip().strip("'").strip('"')

    def _addkey(self, key):
        if self.keylist.count(key) == 0:
            self.keylist.append(key)

    def __setitem__(self, key, value):
        if self.fixed_keys and self.keylist.count(key) == 0:
            raise KeyError("Invalid parameter key: " + key)
        self._addkey(key)
        dict.__setitem__(self, key, value)

    def write(self, outfile=None):
        """
        Write current set of parameters to an output par file. If
        outfile is given, then reset self.filename to the absolute
        path of outfile.
        """
        if outfile is not None:
            self.filename = os.path.abspath(outfile)
        output = open(self.filename, 'w')
        for key in self.keylist:
            try:
                output.write('%s = %s\n' % (key, self[key]))
            except TypeError:
                output.write('%s = %s\n' % (key, repr(self[key])))
        output.close()

    def update(self, pars):
        """
        Reimplement from dict base class to check self.fixed_keys and
        to update self.keylist.
        """
        if not self.fixed_keys:
            for key in pars.keylist:
                self._addkey(key)
                self[key] = pars[key]

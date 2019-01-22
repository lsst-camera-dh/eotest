import os
import subprocess
import lsst.utils


def getVersion():
    import lsst.eotest
    try:
        package_dir = os.path.split(lsst.eotest.__file__)[0]
        command = 'cd %s; git describe --tags' % package_dir
        tag = subprocess.check_output(command, shell=True,
                                      stderr=subprocess.STDOUT)
        return tag.strip()
    except subprocess.CalledProcessError:
        # Not in a git repository, so assume this is an installation
        # from a versioned tarball.
        for item in lsst.eotest.__file__.split(os.path.sep):
            if item.startswith('eotest-'):
                return item[len('eotest-'):]
    return 'unknown'

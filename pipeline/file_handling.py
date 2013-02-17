import os
import subprocess

def setVariable(key, value):
    subprocess.call('pipelineSet %s %s' % (key, value), shell=True)

def get_file_list(prefix):
    numfiles = int(os.environ["NUM%sFILES" % prefix])
    my_files = []
    for i in range(numfiles):
        my_files.append(os.environ["%s_%02i" % (prefix, i)])
    return my_files

def export_file_list(files, prefix):
    setVariable("NUM%sFILES" % prefix, "%s" % len(files))
    for i, item in enumerate(files):
        setVariable("%s_%02i" % (prefix, i), item)

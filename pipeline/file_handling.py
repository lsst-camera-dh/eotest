import os

def get_file_list(prefix):
    numfiles = int(os.environ["NUM%sFILES" % prefix])
    my_files = []
    for i in range(numfiles):
        my_files.append(os.environ["%s_%02i" % (prefix, i)])
    return my_files

def export_file_list(files, prefix):
    import pipeline
    pipeline.setVariable("NUM%sFILES" % prefix, "%s" % len(files))
    for i, item in enumerate(files):
        pipeline.setVariable("%s_%02i" % (prefix, i), item)


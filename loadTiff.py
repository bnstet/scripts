import tifffile as tif
import numpy as np
import os




def tiffToArray(fileOrDir, concatAxis=0):
    if os.path.isdir(fileOrDir):
        tfiles = glob.glob( os.path.join(fileOrDir, '*.tiff')) + glob.glob( os.path.join(fileOrDir, '*.tif'))
        tfiles = sorted(list(set(tfiles)))
        return tif.TiffSequence(tfiles).asarray()
    elif os.path.splitext(fileOrDir)[1] in ['.tif', '.tiff']:
        return tif.imread(fileOrDir)


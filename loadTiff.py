import tifffile as tif
import numpy as np
import os
import glob



def tiffToArray(fileOrDir, concatAxis=0):
    if os.path.isdir(fileOrDir):
        tfiles = glob.glob( os.path.join(fileOrDir, '*.tiff')) + glob.glob( os.path.join(fileOrDir, '*.tif'))
        tfiles = sorted(list(set(tfiles)))
        framelist = []
        for tfile in tfiles:
            frame = tif.imread(tfile)
            shape = frame.shape
            frame = np.expand_dims(frame.reshape((-1,shape[-2],shape[-1]), order='C')[0],axis=0)
            framelist = framelist + [frame]
        return np.concatenate(framelist,axis=0)
    elif os.path.splitext(fileOrDir)[1] in ['.tif', '.tiff']:
        return tif.imread(fileOrDir)


import numpy as np
import tifffile as tif
import glob
import os
import sys

def tifs_to_int(tif_file_list, out_dir):
    mov_list = []
    for f in tif_file_list:
        mov_list = mov_list + [tif.imread(f)]
    totmax = max([np.abs(mov).max() for mov in mov_list])
    for mov,f in zip(mov_list,tif_file_list):
        tif.imsave(data=np.int16(mov*((2**15)/(totmax+10))), file=os.path.join(out_dir,os.path.basename(f)))


if __name__=="__main__":
    tifs_to_int(sys.argv[1], sys.argv[2])

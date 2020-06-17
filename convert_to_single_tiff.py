"""
convert all directories matching "pattern" in "base_dir" to single tif files 

python convert_to_single_tiff.py base_dir pattern

arguments:
pattern: shell-like pattern to match when looking for directories. enclose in quotes
"""

from loadTiff import tiffToArray
import sys,os
import glob
import tifffile as tif



def convert_to_single_tiff(base_dir,pattern):
    dirnames = glob.glob(os.path.join(base_dir,pattern))
    dirnames = [x for x in dirnames if os.path.isdir(x)]
    print('Found {} directories matching pattern {}'.format(len(dirnames),pattern))
    
    for n,dirname in enumerate(dirnames):
        print('({}/{}) converting {} to single .tif'.format(n,len(dirnames),os.path.basename(dirname)))
        mov = tiffToArray(dirname)
        tif.imsave(file=dirname + '.tif',data=mov)

if __name__=="__main__":
    base_dir = sys.argv[1]
    pattern = sys.argv[2]
    convert_to_single_tiff(base_dir,pattern)

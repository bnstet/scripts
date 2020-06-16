"""
run suite2p pipeline

python run_s2p.py ops_file data_path pattern

arguments:
ops_file: path to .npy file containing the "ops" dictionary for suite2p pipeline options
data_path: directory containing the .tif files to be processed
pattern: pattern of tif files to look for. must be sitting in the first level of data_path. enclose in quotes

"""
import sys,os,glob
from suite2p.run_s2p import run_s2p
import numpy as np
from shutil import rmtree

ops_file = sys.argv[1]
data_path = sys.argv[2]

tiff_list = sorted([os.path.basename(x) for x in glob.glob(os.path.join(data_path,pattern))])

ops = np.load(ops_file, allow_pickle=True).item()

db = {
      'h5py': [], # a single h5 file path
      'h5py_key': 'data',
      'look_one_level_down': False, # whether to look in ALL subfolders when searching for tiffs
      'data_path': [data_path], # a list of folders with tiffs 
                                             # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)
      'tiff_list' : tiff_list,
      'subfolders': [], # choose subfolders of 'data_path' to look in (optional)
      'fast_disk': '', # string which specifies where the binary file will be stored (should be an SSD)
    }

if os.path.exists(os.path.join(data_path,'suite2p')):
    rmtree(os.path.join(data_path,'suite2p'))


run_s2p(ops=ops,db=db)

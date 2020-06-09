import sys,os,glob
from suite2p.run_s2p import run_s2p
import numpy as np
from shutil import rmtree

ops_file = sys.argv[1]
data_path = sys.argv[2]

ops = np.load(ops_file, allow_pickle=True).item()

db = {
      'h5py': [], # a single h5 file path
      'h5py_key': 'data',
      'look_one_level_down': False, # whether to look in ALL subfolders when searching for tiffs
      'data_path': [data_path], # a list of folders with tiffs 
                                             # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)
                                            
      'subfolders': [], # choose subfolders of 'data_path' to look in (optional)
      'fast_disk': '', # string which specifies where the binary file will be stored (should be an SSD)
    }

if os.path.exists(os.path.join(data_path,'suite2p')):
    rmtree(os.path.join(data_path,'suite2p'))


run_s2p(ops=ops,db=db)

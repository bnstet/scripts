#!/usr/bin/env python
# -*- coding: utf-8 -*-



import glob
import numpy as np
import os
import logging
import matplotlib.pyplot as plt

try:
    if __IPYTHON__:
        # this is used for debugging purposes only.
        get_ipython().magic('load_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    pass

import caiman as cm
from caiman.paths import caiman_datadir
from caiman.source_extraction import cnmf as cnmf
from caiman.utils.utils import download_demo

import cv2
import logging
import argparse
import tifffile as tif



# movie chunks to be processed
file_chunks = [r'C:\Users\bnste\Documents\scripts\jon_2p_data\JG24831_181210_field1_behavior_00001_00001.tif', 
    r'C:\Users\bnste\Documents\scripts\jon_2p_data\JG24831_181210_field1_behavior_00001_00002.tif',
    r'C:\Users\bnste\Documents\scripts\jon_2p_data\JG24831_181210_field1_behavior_00001_00003.tif']

# channels from scanimage (to be de-interleaved)
n_channels = 2

# which channels to save (suffix will indicate channel) (0-indexed)
saved_channels = [0]

# filter out bad frames with mean frame intensity more than max_dev*(std dev of mean frame intensity) from the average mean frame intensity
# use None to disable filtering
# NOTE: also used to filter base on frame-by-frame differences
max_dev = 4

for n,f in enumerate(file_chunks):
    print('Processing movie file {}'.format(f))
    mov = tif.imread(f)
    if len(mov.shape)==3:
            mov = np.expand_dims(mov,axis=1).astype('float32')
    for ch_num in saved_channels:
        print('Processing channel {}'.format(ch_num))
        ch_out_base, ch_out_file = os.path.split(f)
        pre_ext, ext = os.path.splitext(ch_out_file)
        if max_dev is not None:
            filter_tag = '_filtered'
        else:
            filter_tag = ''
        ch_out_file = os.path.join(ch_out_base , pre_ext) + filter_tag + '_channel_{}'.format(ch_num) + '.tif'
        ch_mov = mov[ch_num::n_channels]
        ch_mean = np.mean(ch_mov)
        ch_means = np.mean(ch_mov,axis=(-1,-2), keepdims=True)
        ch_mean_stddev = np.std(ch_means)


        ch_diff_mean = np.mean(np.diff(ch_mov,axis=0))
        ch_diff_means = np.mean(np.diff(ch_mov,axis=0), axis=(-1,-2), keepdims=True )
        ch_diff_mean_stddev = np.std(ch_diff_means)

        filter_cond = ((ch_means - ch_mean) < max_dev * ch_mean_stddev).flatten() & ( np.concatenate([ [True],  ((ch_diff_means - ch_diff_mean) < max_dev * ch_diff_mean_stddev).flatten() ]))
        
        
        
        ch_mov = ch_mov[filter_cond]
        bad_inds = np.argwhere(~filter_cond)

        print('Removed {} bad frames from file {}, channel {}'.format( len(bad_inds), f, ch_num))
        print('')
        
        if os.path.exists(ch_out_file):
            os.remove(ch_out_file)

        print('Saving to file {}'.format(ch_out_file))
        tif.imsave(ch_out_file, data=ch_mov,  imagej=True)
        print('Save complete.')
        print('')

print('Processing complete.')
        
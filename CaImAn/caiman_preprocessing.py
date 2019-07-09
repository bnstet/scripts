#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Process vids from scanimage during 2p stimulation into files suitable for analysis.

- split out channels into separate vids
- make the vids imagej-compatible
- replace frames with stimulation artifacts with versions interpolated from surrounding good frames

use:
python caiman_preprocessing.py vid_list_file.txt
OR
python caiman_preprocessing.py mov_file.tif

"""


import glob
import numpy as np
import os
import logging
import matplotlib.pyplot as plt
import sys


import caiman as cm

import cv2
import logging
import argparse
import tifffile as tif



# movie chunks to be processed

if '.txt' in sys.argv[1]:
    with open(sys.argv[1], mode='r') as f:
        file_list = [x.replace('\n','') for x in f.readlines()]
else:
    file_list = [sys.argv[1]]


# channels from scanimage (to be de-interleaved)
n_channels = 1

# which channels to save (suffix will indicate channel) (0-indexed)
saved_channels = [0]

# filter out bad frames with a small fraction of extreme values
# the test statistic is LOW for stim frames
# set stim_percentile to be a slight overestimate of the percentage of stim frames in the video
filter_frames = True
stim_percentile = 7

# use an AR filter to smooth out the movie
# alpha is the smoothing parameter (i.e. alpha=1 means a flat moving average of size nSmooth, while alpha=0 means no smoothing)
# use alpha=0 to disable
alpha=0
nSmooth = 1
smoothWindow = np.flip(np.array([alpha**x for x in range(nSmooth)])) + 1e-6
smoothWindow = smoothWindow / smoothWindow.sum()


for n,f in enumerate(file_list):
    print('Processing movie file {}'.format(f), flush=True)
    mov = tif.imread(f)
    if len(mov.shape)==3:
            mov = np.expand_dims(mov,axis=1).astype('float32')
    for ch_num in saved_channels:
        print('Processing channel {}'.format(ch_num), flush=True)
        ch_out_base, ch_out_file = os.path.split(f)
        pre_ext, ext = os.path.splitext(ch_out_file)
        if filter_frames:
            filter_tag = '_filtered'
        else:
            filter_tag = ''
        ch_out_file = os.path.join(ch_out_base , pre_ext) + filter_tag + '_channel_{}'.format(ch_num) + '.tif'
        ch_mov = mov[ch_num::n_channels]
       
        if filter_frames:
            ch_mean = np.mean(ch_mov)
            ch_means = np.mean(ch_mov,axis=(-1,-2), keepdims=True)
            
            frameshape = ch_mov.shape[-2:]
            totpix = frameshape[0]*frameshape[1]

            # get test statistic to measure extreme values: ratio of differences between percentiles
            # will be LOW for stim frames, which have a small fraction of very high values
            pcts = np.percentile(ch_mov, (50, 75, 99.9),axis=(-1,-2))
            mystat =  (pcts[1] - pcts[0]  ) / (pcts[2] - pcts[1])
            mystat = np.nan_to_num(mystat) # nan values sent to 0
            statthresh = np.percentile(mystat,stim_percentile) 
            filter_cond = mystat > statthresh
            
            
            good_inds = np.argwhere(filter_cond)[:,0]
            bad_inds = np.argwhere(~filter_cond)[:,0]
            
            if len(good_inds)==0:
                raise Exception('all frames have been filtered out from file {}'.format(f))
                
            for bad_ind in bad_inds:
                try:
                    prev_good_ind = good_inds[np.argwhere( good_inds < bad_ind)[-1,0]]
                except:
                    prev_good_ind = -np.Inf
                try:
                    next_good_ind = good_inds[np.argwhere( good_inds > bad_ind)[0,0]]
                except:
                    next_good_ind = np.Inf
                if prev_good_ind < 0:
                    ch_mov[bad_ind] = ch_mov[next_good_ind]
                elif next_good_ind > good_inds.max():
                    ch_mov[bad_ind] = ch_mov[prev_good_ind]
                else:
                    gapfrac = (bad_ind - prev_good_ind) / (next_good_ind - prev_good_ind)
                    ch_mov[bad_ind] = ( 1 - gapfrac) * ch_mov[prev_good_ind] + gapfrac * ch_mov[next_good_ind]
                    
                    
                    
                    
            print('Replaced {} bad frames from file {}, channel {}'.format( len(bad_inds), f, ch_num), flush=True)
            print('')

        
        # perform smoothing
        if alpha==0:
            smoothMov=ch_mov
        else:
            smoothMov = np.zeros_like(ch_mov)
            for n, mult in enumerate(smoothWindow):
                movFinalFrame = -(nSmooth-n-1) if nSmooth > n+1 else None
                smoothMov[(nSmooth-n-1):] = smoothMov[(nSmooth-n-1):] + mult * ch_mov[:movFinalFrame]
            initSmoothComp = 1/np.cumsum(np.flip(smoothWindow))

            print('{} {}'.format(smoothMov.shape, initSmoothComp.shape),flush=True)
            smoothMov[:nSmooth] = smoothMov[:nSmooth] * initSmoothComp.reshape((-1,) + tuple(np.ones(len(smoothMov.shape)-1, dtype='int')))
        
        if os.path.exists(ch_out_file):
            os.remove(ch_out_file)

        print('Saving to file {}'.format(ch_out_file), flush=True)
        tif.imsave(ch_out_file, data=smoothMov,  imagej=True)
        print('Save complete.')
        print('', flush=True)

print('Processing complete.', flush=True)
        

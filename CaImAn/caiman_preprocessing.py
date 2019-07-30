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

from PIL import Image


SAVE_TYPE = 'float32'


# handle input arguments
parser = argparse.ArgumentParser(description='filter, smooth, and resize 2p video recordings')
parser.add_argument('file_list', type=str, nargs='+', help='list of files. if one .txt file is provided, that file will be read to retrieve the contained file list (line-separated)')
parser.add_argument('--out_res', type=int, nargs='+', default=None, help='output resolution (defaults to preserving resolution)')
parser.add_argument('--n_channels', type=int, default=1, help='number of input channels')
parser.add_argument('--saved_channels', type=int, nargs='+', default=[0], help='input channels to save (0-indexed)')
parser.add_argument('--pct_thresh', type=float, default=99.99, help='pixel percentile threshold; use 100 for no percentile filtering')
parser.add_argument('--alpha', type=float, default=0, help='AR filter smoothing parameter, alpha=1 means a flat moving average of size n_smooth, while alpha=0 means no smoothing. default 0')
parser.add_argument('--n_smooth', type=int, default=1, help='size of smoothing window. 1 frame is the default')
parser.add_argument('--keep_bad_frames', action='store_true', help='if this flag is activated, keep outlier frames instead of interpolating through them')
parser.add_argument('--ds_factor', type=int, default=1, help='downsample factor; one frame is kept every ds_factor frames. applied prior to smoothing!')
parser.add_argument('--med_ds', dest='ds_func', action='store_const', const=np.median, default=np.mean, help='use median instead of mean for downsampling')


args = parser.parse_args()

print('Using arguments: {}'.format(args),flush=True)

out_res = tuple(args.out_res) if args.out_res is not None else None
n_channels = args.n_channels
saved_channels = args.saved_channels
filter_frames = ~args.keep_bad_frames
pct_thresh = args.pct_thresh
file_list = args.file_list
alpha = args.alpha
nSmooth = args.n_smooth
ds_factor = args.ds_factor
ds_func = args.ds_func

# get file list from text file if text file is in the file_list input variable
if '.txt' in file_list[0]:
    with open(file_list[0], mode='r') as f:
        file_list = [x.replace('\n','') for x in f.readlines()]



# create smoothing window template from parameters
smoothWindow = np.flip(np.array([ (alpha**x if x > 0 else 1) for x in range(nSmooth)])) + 1e-6
smoothWindow = smoothWindow / smoothWindow.sum()



def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh





for n,f in enumerate(file_list):
    print('Processing movie file {}'.format(f), flush=True)
    mov = tif.imread(f).astype(SAVE_TYPE).squeeze()
    for ch_num in saved_channels:
        print('Processing channel {}'.format(ch_num), flush=True)
        ch_out_base, ch_out_file = os.path.split(f)
        pre_ext, ext = os.path.splitext(ch_out_file)
        if filter_frames:
            filter_tag = '_filt'
        else:
            filter_tag = ''
        if out_res is not None:
            res_tag = '_{}x{}'.format(out_res[0], out_res[1])
        else:
            res_tag = '_{}x{}'.format(mov.shape[-2], mov.shape[-1])
        if ds_factor != 1:
            ds_tag = '_ds{}x'.format(ds_factor)
        else:
            ds_tag = ''
        ch_out_file = os.path.join(ch_out_base , pre_ext) + filter_tag + '_chan_{}'.format(ch_num) + res_tag + ds_tag + '.tif'
        ch_mov = mov[ch_num::n_channels]
       

        # resize, if necessary
        if not (ch_mov.shape[:-2] == out_res  or  out_res is None ):
            ch_mov_tmp = np.empty((len(ch_mov),) + out_res)
            for n,frame in enumerate(ch_mov):
                im = Image.fromarray(frame)
                ch_mov_tmp[n] = im.resize(out_res, resample=Image.BICUBIC)
            ch_mov = ch_mov_tmp
            del ch_mov_tmp




        # bound pixel values from above. this will be used to clip the pixels after finding outliers.
        pix_pct_cutoff = np.percentile(ch_mov, pct_thresh)
       

        if filter_frames:
           
            print('Filtering out bad frames...', flush=True)

            frameshape = ch_mov.shape[-2:]
            totpix = frameshape[0]*frameshape[1]

            ch_mov_diff = np.diff(ch_mov,axis=0)

            mystat = np.abs(ch_mov_diff.mean(axis=(-1,-2))) # criterion based on per-frame average of frame changes
            is_bad_ind = is_outlier(mystat)
            is_bad_ind = np.concatenate([is_bad_ind[:-1] | is_bad_ind[1:],[True]])
            
            good_inds = np.argwhere(~is_bad_ind)[:,0]
            bad_inds = np.argwhere(is_bad_ind)[:,0]
            
            if len(good_inds)==0:
                raise Exception('all frames have been filtered out from file {}'.format(f))
           
            # use upper threshold to clip before interpolating through bad frames
            ch_mov = np.clip(ch_mov, a_min=None, a_max=pix_pct_cutoff)

     
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

        # perform downsampling
        print('Downsampling...', flush=True)
        if len(ch_mov) % ds_factor != 0:
            raise Exception('Movie length {} is not divisible by ds_factor {}'.format(len(ch_mov), ds_factor))
        
        ch_mov = ch_mov.reshape((len(ch_mov) // ds_factor, ds_factor) + ch_mov.shape[1:], order='C')
        ch_mov = ds_func(ch_mov, axis=1)
        
        
        # perform smoothing
        print('Smoothing...', flush=True)
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
        tif.imsave(ch_out_file, data=np.expand_dims(smoothMov.astype(SAVE_TYPE), 1),  imagej=True)
        print('Save complete.')
        print('', flush=True)

print('Processing complete.', flush=True)
        

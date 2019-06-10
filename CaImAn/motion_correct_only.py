#!/usr/bin/env python
"""
Motion correction via CaImAn.

copyright GNU General Public License v2.0
authors: @agiovann and @epnev
"""

print('Loading libraries')

import cv2
import glob
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import tifffile as tif

try:
    cv2.setNumThreads(0)
except:
    pass

try:
    if __IPYTHON__:
        # this is used for debugging purposes only. allows to reload classes
        # when changed
        get_ipython().magic('load_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    pass


import caiman as cm
from caiman.motion_correction import MotionCorrect
from caiman.source_extraction.cnmf import params as params
from caiman.utils.utils import download_demo



parser = argparse.ArgumentParser(description='Execute the CaImAn (Calcium Imaging Analysis) motion correction operation on a video recording (.tif file or .tiff stack).')
parser.add_argument('infile', type=str, nargs=1, help='file name (for .tif file) or folder (for .tiff stack)')
parser.add_argument('--outfile',type=str, nargs=1, default=[], help='file name under which to store the output motion correction movie')
parser.add_argument('--mc_temp',type=str, nargs='?', default=None, help='template .tif file for motion correction')
args = vars(parser.parse_args())

infile = args['infile']
outfile = args['outfile']
mc_temp = args['mc_temp']

print('Using arguments list: {}'.format(args))

orig_in = infile[0]


if '.txt' in infile[0]:
    with open(infile[0], 'r') as f:
        flist = []
        for line in f:
            flist = flist + [line.replace('\n','')]
    fname = flist
elif os.path.isdir(infile[0]):
    tmpMovPath = os.path.join(infile[0], 'tmp_mov.tif')
    if not os.path.isfile(tmpMovPath):
        tfiles = glob.glob( os.path.join(infile[0], '*.tiff')) + glob.glob( os.path.join(infile[0], '*.tif'))
        tfiles = sorted(list(set(tfiles)))
        with tif.TiffWriter(tmpMovPath, bigtiff=True) as writer:
            for i in range(len(tfiles)):
                if (i%100) == 0:
                    print("Writing movie to temp file, frame {}/{}".format(i,len(tfiles)))
                im = tif.imread(tfiles[i]).astype('float16')
                if len(im.shape)==2:
                    writer.save(im , contiguous=True, compress=6, photometric='minisblack')
                else:
                    for j in range(len(im)):
                        writer.save(im[j], contiguous=True, compress=6, photometric='minisblack')
    fname = [tmpMovPath]
else:
    tmpMovPath = None
    fname = infile

# set outfile name
outfileStore = outfile
outfile = []
in_file_ext = []
for orig_in in fname:
    if len(outfileStore)==0:
        if os.path.isdir(orig_in):
            outfile = outfile + [orig_in + '_mc']
            in_file_ext = in_file_ext + ['']
        else:
            in_filename, in_file_ext_new = os.path.splitext(orig_in)
            in_file_ext = in_file_ext + [in_file_ext_new]
            outfile = outfile + [in_filename + '_mc' + in_file_ext_new]
    else:
        outfile = outfileStore


# load in template file
if mc_temp is not None:
    print('Loading motion correction reference template {}'.format(mc_temp),flush=True)
    mc_temp = tif.imread(mc_temp)
    if len(mc_temp.shape)>2:
        mc_temp = mc_temp[0]


#%%

def main():
    pass  # For compatibility between running under Spyder and the CLI

    #%% Select file(s) to be processed (download if not present)
    fnames = fname

    #%% First setup some parameters for data and motion correction
    
    n_processes = 18

    # dataset dependent parameters

    nChannels = 1
    channelToKeep = 0

    fr = 30            # imaging rate in frames per second
    decay_time = 0.4    # length of a typical transient in seconds
    dxy = (2., 2.)      # spatial resolution in x and y in (um per pixel)
    # note the lower than usual spatial resolution here
    max_shift_um = (12., 12.)       # maximum shift in um
    patch_motion_um = (100., 100.)  # patch size for non-rigid correction in um

    # motion correction parameters
    pw_rigid = False       # flag to select rigid vs pw_rigid motion correction
    # maximum allowed rigid shift in pixels
    #max_shifts = [int(a/b) for a, b in zip(max_shift_um, dxy)]
    max_shifts = (20,20)
    # start a new patch for pw-rigid motion correction every x pixels
    #strides = tuple([int(a/b) for a, b in zip(patch_motion_um, dxy)])
    strides = (32,32)
    # overlap between pathes (size of patch in pixels: strides+overlaps)
    overlaps = (18, 18)
    # maximum deviation allowed for patch with respect to rigid shifts
    max_deviation_rigid = 5


  # estimate minimum of movie
    min_mov = np.min(tif.imread(fname[0]))

    mc_dict = {
        'fnames': fnames,
        'fr': fr,
        'decay_time': decay_time,
        'dxy': dxy,
        'pw_rigid': pw_rigid,
        'max_shifts': max_shifts,
        'strides': strides,
        'overlaps': overlaps,
        'max_deviation_rigid': max_deviation_rigid,
        'border_nan': 'min',
        'min_mov': min_mov
    }

    opts = params.CNMFParams(params_dict=mc_dict)

    
    # %% start a cluster for parallel processing
    c, dview, n_processes = cm.cluster.setup_cluster(
        backend='local', n_processes=n_processes, single_thread=False)

    print('checkpoint 1: mcorrect start', flush=True)


    # %%% MOTION CORRECTION
    # first we create a motion correction object with the specified parameters
    mc = MotionCorrect(fname, dview=dview, **opts.get_group('motion'))
    # note that the file is not loaded in memory

    # %% Run (piecewise-rigid motion) correction using NoRMCorre
    mc.motion_correct(save_movie=True, template=mc_temp)
    border_to_0 = 0 if mc.border_nan is 'copy' else mc.border_to_0
 
    if pw_rigid:
        mmap_files = mc.fname_tot_els
    else:
        mmap_files = mc.fname_tot_rig

    print('mmap file list: {}'.format(mmap_files), flush=True)

    for n,mmap_file in enumerate(mmap_files):
        fname_new = cm.save_memmap([mmap_file], base_name='memmap_{}_'.format(n), order='C',border_to_0=border_to_0)
    
        (mov, frame_dims, num_frames) = cm.load_memmap(fname_new)
        mov = mov.reshape( list(frame_dims) + [num_frames]).transpose((2,1,0))
        mov = mov[channelToKeep::nChannels]


        # currently format output based on input (stack or single-file)
        #if in_file_ext[n] in ['.tif', '.tiff']:
        if '.tif' in outfile[n]:
            tif.imsave(outfile[n], mov, imagej=True)
        else:
            print("saving output {}".format(outfile[n]), flush=True)
            if not os.path.exists(outfile[n]):
                os.mkdir(outfile[n])
            for i in range(len(mov)):
                frame = mov[i]
                frame_file = os.path.join(outfile[n], 'frame{:06d}.tiff'.format(i))
                tif.imsave(frame_file, frame, imagej=True)
    
    print('Save complete. Removing temporary files.', flush=True)

    if tmpMovPath is not None:
        os.remove(tmpMovPath)
    [os.remove(x) for x in mmap_files];


    print('Motion correction script complete.', flush=True)

    


# %%
# This is to mask the differences between running this demo in Spyder
# versus from the CLI
if __name__ == "__main__":
    main()

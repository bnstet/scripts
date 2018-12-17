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
parser.add_argument('--outfile',type=str, nargs=1, default=[], help='file name under which to store the output hdf5 file (use .hdf5 suffix)')
args = vars(parser.parse_args())

infile = args['infile']
outfile = args['outfile']

print('Using arguments list: {}'.format(args))

orig_in = infile[0]

if os.path.isdir(infile[0]):
    tmpMovPath = os.path.join(infile[0], 'tmp_mov.tif')
    if not os.path.isfile(tmpMovPath):
        tfiles = glob.glob( os.path.join(infile[0], '*.tiff')) + glob.glob( os.path.join(infile[0], '*.tif'))
        tfiles = sorted(list(set(tfiles)))
        with tif.TiffWriter(tmpMovPath, bigtiff=True) as writer:
            for i in range(len(tfiles)):
                if (i%100) == 0:
                    print("Writing movie to temp file, frame {}/{}".format(i,len(tfiles)))
                im = tif.imread(tfiles[i])
                if len(im.shape)==2:
                    writer.save(im, compress=6, photometric='minisblack')
                else:
                    for j in range(len(im)):
                        writer.save(im[j], compress=6, photometric='minisblack')
    fname = [tmpMovPath]
else:
    fname = infile

in_filename, in_file_ext = os.path.splitext(orig_in)
if len(outfile)==0:
    outfile = in_filename + '_mc' + in_file_ext
else:
    outfile = outfile[0]

#%%

def main():
    pass  # For compatibility between running under Spyder and the CLI

    #%% Select file(s) to be processed (download if not present)
    fnames = fname

    #%% First setup some parameters for data and motion correction
    
    n_processes = 12

    # dataset dependent parameters
    fr = 3.6             # imaging rate in frames per second
    decay_time = 0.4    # length of a typical transient in seconds
    dxy = (2., 2.)      # spatial resolution in x and y in (um per pixel)
    # note the lower than usual spatial resolution here
    max_shift_um = (12., 12.)       # maximum shift in um
    patch_motion_um = (100., 100.)  # patch size for non-rigid correction in um

    # motion correction parameters
    pw_rigid = False       # flag to select rigid vs pw_rigid motion correction
    # maximum allowed rigid shift in pixels
    max_shifts = [int(a/b) for a, b in zip(max_shift_um, dxy)]
    # start a new patch for pw-rigid motion correction every x pixels
    strides = tuple([int(a/b) for a, b in zip(patch_motion_um, dxy)])
    # overlap between pathes (size of patch in pixels: strides+overlaps)
    overlaps = (24, 24)
    # maximum deviation allowed for patch with respect to rigid shifts
    max_deviation_rigid = 3

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
        'border_nan': 'copy'
    }

    opts = params.CNMFParams(params_dict=mc_dict)

    
    # %% start a cluster for parallel processing
    c, dview, n_processes = cm.cluster.setup_cluster(
        backend='local', n_processes=n_processes, single_thread=False)

    print('checkpoint 1: mcorrect start')

    # %%% MOTION CORRECTION
    # first we create a motion correction object with the specified parameters
    mc = MotionCorrect(fname[0], dview=dview, **opts.get_group('motion'))
    # note that the file is not loaded in memory

    # %% Run (piecewise-rigid motion) correction using NoRMCorre
    mc.motion_correct(save_movie=True)

    print('checkpoint 2: mcorrect finish; temp file stored in {}'.format(mc.mmap_file[0]))
    print('Save corrected tif stack to {}'.format(outfile))
    
    (mov, frame_dims, num_frames) = cm.load_memmap(mc.mmap_file[0])
    mov = mov.reshape( frame_dims + (num_frames,), order='C')
    # mov is currently in (height, width, nframes format)


    # currently format output based on input (stack or single-file)
    if in_file_ext in ['.tif', '.tiff']:
        tif.imsave(outfile, mov, imagej=True)
    else:
        if not os.path.exists(outfile):
            os.mkdir(outfile)
        mov = mov.transpose((2,1,0)) #mov now has nframes in the first dimension
        for i in range(len(mov)):
            frame = mov[i]
            frame_file = os.path.join(outfile, 'frame{:06d}.tiff'.format(i))
            tif.imsave(frame_file, frame, imagej=True)
    
    print('Save complete. Removing temporary files.')

    os.remove(tmpMovPath)
    os.remove(mc.mmap_file[0])


    print('Motion correction script complete.')

    exit()
    
    # %% MEMORY MAPPING
    border_to_0 = 0 if mc.border_nan is 'copy' else mc.border_to_0
    # you can include the boundaries of the FOV if you used the 'copy' option
    # during motion correction, although be careful about the components near
    # the boundaries

    # memory map the file in order 'C'
    fname_new = cm.save_memmap(mc.mmap_file, base_name='memmap_', order='C',
                                border_to_0=border_to_0)  # exclude borders

    # now load the file
    Yr, dims, T = cm.load_memmap(fname_new)
    images = np.reshape(Yr.T, [T] + list(dims), order='F')
    # load frames in python format (T x X x Y)

    # %% restart cluster to clean up memory
    cm.stop_server(dview=dview)
   

    


# %%
# This is to mask the differences between running this demo in Spyder
# versus from the CLI
if __name__ == "__main__":
    main()

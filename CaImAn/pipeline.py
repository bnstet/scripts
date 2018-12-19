#!/usr/bin/env python
"""
Complete demo pipeline for processing two photon calcium imaging data using the
CaImAn batch algorithm. The processing pipeline included motion correction,
source extraction and deconvolution. The demo shows how to construct the
params, MotionCorrect and cnmf objects and call the relevant functions. You
can also run a large part of the pipeline with a single method (cnmf.fit_file)
See inside for details.

Demo is also available as a jupyter notebook (see demo_pipeline.ipynb)
Dataset couresy of Sue Ann Koay and David Tank (Princeton University)

This demo pertains to two photon data. For a complete analysis pipeline for
one photon microendoscopic data see demo_pipeline_cnmfE.py

copyright GNU General Public License v2.0
authors: @agiovann and @epnev
"""

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
from caiman.source_extraction.cnmf import cnmf as cnmf
from caiman.source_extraction.cnmf import params as params
from caiman.utils.utils import download_demo



parser = argparse.ArgumentParser(description='Execute the CaImAn (Calcium Imaging Analysis) pipeline on a video recording (.tif file or .tiff stack). Motion correct, determine neuronal regions, extract signals and denoise.')
parser.add_argument('infile', type=str, nargs=1, help='file name (for .tif file) or folder (for .tiff stack)')
parser.add_argument('outfile',type=str, nargs=1, help='file name under which to store the output hdf5 file (use .hdf5 suffix)')
parser.add_argument('--config', type=str, nargs='?', help='file name for config file. should consist of a definition of a list of python dictionaries' )
parser.add_argument('--log_fname', type=str, nargs='?', default="logs/caiman_processing.log", help='file name under which to save a progress log. leave out to save to caiman_processing.log')
parser.add_argument('--mc_fname', type=str, nargs='?', default="", help='file name under which to save motion-corrected video. leave out to not save (default behavior)')
parser.add_argument('--nomc', action='store_true', help='if used, then no motion correction will be run on the input video')
parser.add_argument('--slurmid', type=int, nargs='?', default=0, help='slurm ID for multi-node processing. leave out to default to 0')
args = vars(parser.parse_args())

infile = args['infile']
outfile = args['outfile'][0]
log_fname = args['log_fname']
mc_fname = args['mc_fname']
nomc = args['nomc']
slurmid = args['slurmid']

splitlog = log_fname.split('.')
head = str.join('.',splitlog[:-1])
tail = '.' + splitlog[-1]
log_fname = head + '_' + str(slurmid) + tail


# %%
# Set up the logger (optional); change this if you like.
# You can log to a file using the filename parameter, or make the output more
# or less verbose by setting level to logging.DEBUG, logging.INFO,
# logging.WARNING, or logging.ERROR

logging.basicConfig(format=
                    "%(relativeCreated)12d [%(filename)s:%(funcName)20s():%(lineno)s]"\
                    "[%(process)d] %(message)s",
                    level=logging.INFO, filename = log_fname)

#%%

# save single-file movie if given folder of tiffs
"""
if os.path.isdir(infile[0]):
    tmpMovPath = os.path.join(infile[0], 'tmp_mov.hdf5')
    cm.load(glob.glob( os.path.join(infile[0], '*.tiff'))).save(tmpMovPath)
    fname = [tmpMovPath]
else:
    fname = infile
"""



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
                im = tif.imread(tfiles[i])
                if len(im.shape)==2:
                    writer.save(im, compress=6, photometric='minisblack')
                else:
                    for j in range(len(im)):
                        writer.save(im[j], compress=6, photometric='minisblack')
    fname = [tmpMovPath]
else:
    fname = infile

print('Using file(s): {}'.format(fname))

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

    # %% play the movie (optional)
    # playing the movie using opencv. It requires loading the movie in memory.
    # To close the video press q
    display_images = False

    if display_images:
        m_orig = cm.load_movie_chain(fnames)
        ds_ratio = 0.2
        moviehandle = m_orig.resize(1, 1, ds_ratio)
        moviehandle.play(q_max=99.5, fr=60, magnification=2)

    # %% start a cluster for parallel processing
    c, dview, n_processes = cm.cluster.setup_cluster(
        backend='local', n_processes=n_processes, single_thread=False)

    print('checkpoint 1: mcorrect')

    # %%% MOTION CORRECTION
    # first we create a motion correction object with the specified parameters
    mc = MotionCorrect(fnames, dview=dview, **opts.get_group('motion'))
    # note that the file is not loaded in memory

    # %% Run (piecewise-rigid motion) correction using NoRMCorre
    mc.motion_correct(save_movie=True)

    # %% compare with original movie
    if display_images:
        m_orig = cm.load_movie_chain(fnames)
        m_els = cm.load(mc.mmap_file)
        ds_ratio = 0.2
        moviehandle = cm.concatenate([m_orig.resize(1, 1, ds_ratio) - mc.min_mov*mc.nonneg_movie,
                                        m_els.resize(1, 1, ds_ratio)], axis=2)
        moviehandle.play(fr=60, q_max=99.5, magnification=2)  # press q to exit

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
    c, dview, n_processes = cm.cluster.setup_cluster(
        backend='local', n_processes=None, single_thread=False)
    n_processes = int(n_processes)

    
    # %%  parameters for source extraction and deconvolution
    p = 1                    # order of the autoregressive system
    gnb = 2                  # number of global background components
    merge_thresh = 0.8       # merging threshold, max correlation allowed
    rf = 15                  # half-size of the patches in pixels. e.g., if rf=25, patches are 50x50
    stride_cnmf = 4          # amount of overlap between the patches in pixels
    K = 4                    # number of components per patch
    gSig = [4, 4]            # expected half size of neurons in pixels
    method_init = 'greedy_roi'   # initialization method (if analyzing dendritic data using 'sparse_nmf')
    ssub = 2                     # spatial subsampling during initialization
    tsub = 2                     # temporal subsampling during intialization

    # parameters for component evaluation
    opts_dict = {'fnames': fnames,
                    'fr': fr,
                    'nb': gnb,
                    'rf': rf,
                    'K': K,
                    'gSig': gSig,
                    'stride': stride_cnmf,
                    'method_init': method_init,
                    'rolling_sum': True,
                    'merge_thr': merge_thresh,
                    'n_processes': n_processes,
                    'only_init': True,
                    'ssub': ssub,
                    'tsub': tsub}

    opts.change_params(params_dict=opts_dict)



    print('checkpoint 2: patch cnmf')

    # %% RUN CNMF ON PATCHES
    # First extract spatial and temporal components on patches and combine them
    # for this step deconvolution is turned off (p=0)

    opts.set('temporal', {'p': 0})
    cnm = cnmf.CNMF(n_processes, params=opts, dview=dview)
    cnm = cnm.fit(images)

    # %% ALTERNATE WAY TO RUN THE PIPELINE AT ONCE
    #   you can also perform the motion correction plus cnmf fitting steps
    #   simultaneously after defining your parameters object using
    #  cnm1 = cnmf.CNMF(n_processes, params=opts, dview=dview)
    #  cnm1.fit_file(motion_correct=True)

    # %% plot contours of found components
    Cn = cm.local_correlations(images, swap_dim=False)
    Cn[np.isnan(Cn)] = 0
    #cnm.estimates.plot_contours(img=Cn)
    #plt.title('Contour plots of found components')


    print('checkpoint 3: eval components')

    # %% RE-RUN seeded CNMF on accepted patches to refine and perform deconvolution
    cnm.params.set('temporal', {'p': p})
    cnm2 = cnm.refit(images, dview=dview)
    
    # %% COMPONENT EVALUATION
    # the components are evaluated in three ways:
    #   a) the shape of each component must be correlated with the data
    #   b) a minimum peak SNR is required over the length of a transient
    #   c) each shape passes a CNN based classifier
    min_SNR = 1 # signal to noise ratio for accepting a component
    rval_thr = 0.7  # space correlation threshold for accepting a component
    cnn_thr = 0.97  # threshold for CNN based classifier
    cnn_lowest = 0.1 # neurons with cnn probability lower than this value are rejected

    cnm2.params.set('quality', {'decay_time': decay_time,
                                'min_SNR': min_SNR,
                                'rval_thr': rval_thr,
                                'use_cnn': True,
                                'min_cnn_thr': cnn_thr,
                                'cnn_lowest': cnn_lowest})
    cnm2.estimates.evaluate_components(images, cnm2.params, dview=dview)
    
    # %% PLOT COMPONENTS
    #cnm2.estimates.plot_contours(img=Cn, idx=cnm2.estimates.idx_components)

    
    # %% VIEW TRACES (accepted and rejected)

    if display_images:
        cnm2.estimates.view_components(images, img=Cn,
                                        idx=cnm2.estimates.idx_components)
        cnm2.estimates.view_components(images, img=Cn,
                                        idx=cnm2.estimates.idx_components_bad)
    
    #%% update object with selected components
    cnm2.estimates.select_components(use_object=True)
    
    #%% Extract DF/F values
    cnm2.estimates.detrend_df_f(quantileMin=8, frames_window=250)

    #%% Show final traces
    #cnm2.estimates.view_components(img=Cn)

    #%% reconstruct denoised movie (press q to exit)
    if display_images:
        cnm2.estimates.play_movie(images, q_max=99.9, gain_res=2,
                                    magnification=2,
                                    bpx=border_to_0,
                                    include_bck=False)  # background not shown

    print('checkpoint 4: save')

    #%% save results
    cnm2.save(outfile)

    #%% STOP CLUSTER and clean up log files
    cm.stop_server(dview=dview)
    #log_files = glob.glob('*_LOG_*')
    #for log_file in log_files:
    #    os.remove(log_file)

    print('checkpoint 5: final')



# %%
# This is to mask the differences between running this demo in Spyder
# versus from the CLI
if __name__ == "__main__":
    main()

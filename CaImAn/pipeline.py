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
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import tifffile as tif
import pickle

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

import yaml

from scipy.io import loadmat


parser = argparse.ArgumentParser(description='Execute the CaImAn (Calcium Imaging Analysis) pipeline on a video recording (.tif file or .tiff stack). Motion correct, determine neuronal regions, extract signals and denoise.')

parser.add_argument('infile', type=str, nargs=1, help='file name (for .tif file) or folder (for .tiff stack). can also be a .txt file containing a list of .tif files for multiple chunks.')
parser.add_argument('outfile',type=str, nargs=1, help='file name under which to store the output hdf5 file (use .hdf5 suffix)')
parser.add_argument('--mc_temp',type=str, nargs='?', default=None, help='template .tif file for motion correction')
parser.add_argument('--config', type=str, nargs='?', help='file name for config file. should consist of a definition of a list of python dictionaries' )
parser.add_argument('--log_fname', type=str, nargs='?', default="logs/caiman_processing.log", help='file name under which to save a progress log. leave out to save to caiman_processing.log')
parser.add_argument('--mc_fname', type=str, nargs='?', default="", help='file name under which to save motion-corrected video. leave out to not save (default behavior)')
parser.add_argument('--nomc', action='store_true', help='if used, then no motion correction will be run on the input video/mmap file')
parser.add_argument('--slurmid', type=int, nargs='?', default=0, help='slurm ID for multi-node processing. leave out to default to 0')
parser.add_argument('--Ainfile',type=str, nargs=1, default=None, help='optional file for input masks. should either be .mat or .npy file containing a(n_masks)x(dim_1)x(dim_2) array of masks called "Ain"' )
parser.add_argument('--settingsFile', nargs=1, default='settings/caiman_settings_default.yml', help='YAML file containing CaImAn model parameters')

args = vars(parser.parse_args())

infile = args['infile']
outfile = args['outfile'][0]
mc_temp = args['mc_temp']
log_fname = args['log_fname']
mc_fname = args['mc_fname']
nomc = args['nomc']
slurmid = args['slurmid']

# handle settings file
settingsFile = args['settingsFile']
if type(settingsFile) is list:
    settingsFile = settingsFile[0]
if not os.path.isfile(settingsFile):
    raise FileNotFoundError('settingsFile {} not found!'.format(settingsFile))
try:
    settingsFileDict = yaml.safe_load(open(settingsFile, 'r'))
except yaml.YAMLError as exc:
    print('Error parsing YAML settingsFile {} : {}'.format(settingsFile,exc))
    raise exc


# handle Afile (input mask file)
Ainfile = args['Ainfile']
if Ainfile is not None:
    if type(Ainfile) is list:
        Ainfile = Ainfile[0]
    if '.mat' in Ainfile:
        Ahandle = loadmat(Ainfile)
        Ain = np.squeeze(Ahandle['Ain'])
    elif ('.npy' in Ainfile) or ('.npz' in Ainfile):
        Ahandle = np.load(Ainfile)
        Ain = np.squeeze(Ahandle['Ain'])
    else:
        raise Exception('Ainfile argument must point to a .mat, .npy, or .npz file')
else:
    Ain = None


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

# process input file

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
# get file extension (of first file)
_, fext = os.path.splitext(fname[0])

## process template into ndarray

if mc_temp is not None:
    mc_temp = tif.imread(mc_temp)
    if len(mc_temp.shape)>2:
        mc_temp = mc_temp[0]






#%% Select file(s) to be processed (download if not present)
fnames = fname

#%% First setup some parameters for data and motion correction
    
n_processes = 12

# dataset dependent parameters
fr = 30             # imaging rate in frames per second
decay_time = 1.5    # length of a typical transient in seconds
dxy = (1., 1.)      # spatial resolution in x and y in (um per pixel)
# note the lower than usual spatial resolution here
max_shift_um = (12., 12.)       # maximum shift in um
patch_motion_um = (100., 100.)  # patch size for non-rigid correction in um

# motion correction parameters
pw_rigid = True       # flag to select rigid vs pw_rigid motion correction
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
    'border_nan': 'min'
}

#opts = params.CNMFParams(params_dict=mc_dict)
opts = params.CNMFParams(params_dict=settingsFileDict)

    
# %% start a cluster for parallel processing
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=n_processes, single_thread=False)

print('checkpoint 1: mcorrect', flush=True)

# %%% MOTION CORRECTION
# first we create a motion correction object with the specified parameters
if nomc:
    print('Skipping motion correction due to -nomc flag', flush=True)
    print('saving movies as mmap file', flush=True)
    fname_new = cm.save_memmap(fname, base_name='memmap_{}_full_'.format(slurmid), order='C',
                               border_to_0=0)
    
else:
    mc = MotionCorrect(fnames, dview=dview, **opts.get_group('motion'))
    # note that the file is not loaded in memory
    
    # %% Run (piecewise-rigid motion) correction using NoRMCorre
    mc.motion_correct(save_movie=True, template=mc_temp)
    
        
    # %% MEMORY MAPPING
    border_to_0 = 0 if mc.border_nan is 'copy' else mc.border_to_0
    # you can include the boundaries of the FOV if you used the 'copy' option
    # during motion correction, although be careful about the components near
    # the boundaries
    # memory map the file in order 'C'
    
    save_name = mc.fname_tot_els if opts.get('motion','pw_rigid') else mc.fname_tot_rig
    
    print('saving movies as mmap file', flush=True)
    
    fname_new = cm.save_memmap(save_name, base_name='memmap_{}_full_'.format(slurmid), order='C',
                               border_to_0=border_to_0)  # exclude borders
            
            
            
            
            
        


# reload full mmap sequence for processing

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
gnb = 2                  # number of global background components, if positive, otherwise ring model with settings
# gnb=0 : return background as b and W
# gnb=-1 : retyrb full rank background B
# gnb<-1: don't return background
#Ain = None          # possibility to seed with predetermined binary masks
merge_thr = 0.7       # merging threshold, max correlation allowed
rf = 40                  # half-size of the patches in pixels. e.g., if rf=25, patches are 50x50
stride_cnmf = 20          # amount of overlap between the patches in pixels
K = 8                    # upper bound on components per patch
gSig = [5, 5]            # gaussian width of a 2D gaussian kernel, which approximates a neuron
gSiz = [21, 21]     # average diameter of a neuron, in general 4*gSig+1
bord_px=20 # number of pixels to not consider in the borders)

method_init = 'greedy_roi'   # initialization method (if analyzing dendritic data using 'sparse_nmf'), for standard 2p use greedy_roi, for 1p use corr_pnr


ssub = 1                     # spatial subsampling during initialization, increase if you have memory problems
tsub = 2                     # temporal subsampling during intialization, increase if you have memory problems
center_psf=False      # set to true if there is strong background
low_rank_background = True  # None leaves background of each patch intact, True performs global low-rank approximation if gnb>0
nb_patch = 1        # number of background components (rank) per patch if gnb>0, else it is set automatically (default 1)
min_corr = .6       # min peak value from correlation image (for corr_pnr)
min_pnr = 6        # min peak to noise ration from PNR image (for corr_pnr)
ssub_B = 2          # additional downsampling factor in space for background (for corr_pnr)
ring_size_factor = 1.4  # radius of ring is gSiz*ring_size_factor

only_init = True 

update_background_components=True # sometimes setting to False improve the results


opts_dict={'method_init': method_init,  
           'K': K,
           'gSig': gSig,
           'gSiz': gSiz,
           'merge_thr': merge_thr,
           'p': p,
           'tsub': tsub,
           'ssub': ssub,
           'rf': rf,
           'stride': stride_cnmf,
           'only_init': only_init,    # set it to True to run CNMF-E
           'nb': gnb,
           'nb_patch': nb_patch,
           'method_deconvolution': 'oasis',       # could use 'cvxpy' alternatively
           'low_rank_background': low_rank_background,
           'update_background_components': True,  # sometimes setting to False improve the results
           'min_corr': min_corr,
           'min_pnr': min_pnr,
           'normalize_init': False,               # just leave as is
           'center_psf': center_psf,   
           'ssub_B': ssub_B,
           'ring_size_factor': ring_size_factor,
           'del_duplicates': True,                # whether to remove duplicates from initialization
           'border_pix': bord_px}            


# parameters for component evaluation
#opts_dict = {'fnames': fnames,
#             'fr': fr,
#             'nb': gnb,
#             'rf': rf,
#             'K': K,
#             'gSig': gSig,
#             'stride': stride_cnmf,
#             'method_init': method_init,
#             'rolling_sum': True,
#             'merge_thr': merge_thresh,
#             'n_processes': n_processes,
#             'only_init': True,
#             'ssub': ssub,
#             'tsub': tsub,
#             'center_psf':center_psf
#         }

#opts.change_params(params_dict=opts_dict)


print('checkpoint 2: patch cnmf', flush=True)


# % INITIALIZE CNMF COMPONENTS

#cnm = cnmf.CNMF(n_processes, params=opts, dview=dview)
#cnm.params.set('init',{})




# %% RUN CNMF ON PATCHES
# First extract spatial and temporal components on patches and combine them
# for this step deconvolution is turned off (p=0)


if method_init is "greedy_roi": # standard case
    opts.change_params( {'p': 0})
    cnm = cnmf.CNMF(n_processes, params=opts, Ain=Ain, dview=dview)
    cnm = cnm.fit(images)


    # %% RE-RUN seeded CNMF on accepted patches to refine and perform deconvolution
    cnm.params.change_params( {'p': p})
    cnm2 = cnm.refit(images, dview=dview)
else: # cnmf-E case
    cnm2 = cnmf.CNMF(n_processes, params=opts, Ain=Ain, dview=dview)
    cnm2.fit(images)


print('checkpoint 3: eval components', flush=True)


# %% COMPONENT EVALUATION
# the components are evaluated in three ways:
#   a) the shape of each component must be correlated with the data
#   b) a minimum peak SNR is required over the length of a transient
#   c) each shape passes a CNN based classifier
min_SNR = 1 # signal to noise ratio for accepting a component
rval_thr = 0.7  # space correlation threshold for accepting a component; lower -> more components accepted
use_cnn = True # use CNN classifier on spatial components
cnn_thr = 0.99  # threshold for CNN based classifier
cnn_lowest = 0.1 # neurons with cnn probability lower than this value are rejected

#cnm2.params.set('quality', {'decay_time': decay_time,
#                            'min_SNR': min_SNR,
#                            'rval_thr': rval_thr,
#                            'use_cnn': use_cnn,
#                            'min_cnn_thr': cnn_thr,
#                            'cnn_lowest': cnn_lowest})
cnm2.estimates.evaluate_components(images, cnm2.params, dview=dview)

#%% Extract DF/F values
if cnm2.estimates.b is not None:
    cnm2.estimates.detrend_df_f(quantileMin=8, frames_window=250)

    
#%% update object with selected components
#cnm2.estimates.select_components(use_object=True)


    
print('checkpoint 4: save', flush=True)

#%% save results
cnm2.save(outfile)

#%% STOP CLUSTER and clean up log files
cm.stop_server(dview=dview)
#log_files = glob.glob('*_LOG_*')
#for log_file in log_files:
#    os.remove(log_file)

print('checkpoint 5: final', flush=True)


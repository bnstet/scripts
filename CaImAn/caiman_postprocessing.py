import numpy as np
import caiman as cm
import tifffile as tif

from caiman.source_extraction.cnmf import cnmf as cnmf
from caiman.source_extraction.cnmf import params as params
from caiman.utils.visualization import plot_contours

import matplotlib.pyplot as plt

from scipy.io import savemat

import os
import sys
import glob

from PIL import Image

import twop_rec_tools as tp

import argparse



parser = argparse.ArgumentParser(description='Utilize CaImAn hdf5 analysis output to produce summary images and masks')
parser.add_argument('infile', type=str, nargs=1, help='input hdf5 file')
parser.add_argument('outdir',type=str, nargs=1, help='directory in which to store outputs')
parser.add_argument('--online', action='store_true', help='use to indicate that data is coming from the caiman online (ONACID) algorithm')
parser.add_argument('--nsources', type=int, default=20, help='number of neuronal regions for masks (ranked by decreasing signal stddev)')
parser.add_argument('--maskthresh', type=float, default=.75, help='threshold for masks (lower threshold increases the area)')

args = parser.parse_args()

hdf5File = args.infile[0]
outDir = args.outdir[0].replace('\\', '/')
online = args.online
nSources = args.nsources
maskThresh = args.maskthresh


inBasename, _ = os.path.splitext(os.path.basename(hdf5File))
outArrayFile = os.path.join(outDir, inBasename) + '_caiman_arrays.npz'
outMatlabFile = os.path.join(outDir, inBasename) + '_caiman_arrays.mat'
outImageFolder = os.path.join(outDir, inBasename + '_images')


if not os.path.exists(outDir):
    os.mkdir(outDir)
if not os.path.exists(outImageFolder):
    os.mkdir(outImageFolder)


# load in data

if args.online:
    cnm = cnmf.online_cnmf.load_OnlineCNMF(hdf5File)
else:
    cnm = cnmf.load_CNMF(hdf5File)

if len(cnm.estimates.dims)==1:
    cnm.estimates.dims = cnm.estimates.dims*2

A, C, b, f = tp.cellInfoCaimanHdf5(hdf5File)

meanImageSources = np.tensordot(C.mean(axis=0), A, axes=(0,0)).transpose().astype('float32')
meanImage = meanImageSources + np.tensordot(f.mean(axis=0), b, axes=(0,0)).transpose().astype('float32')

stdImageSources = np.tensordot(C.std(axis=0), A, axes=(0,0)).transpose().astype('float32')
stdImage = stdImageSources + np.tensordot(f.std(axis=0), b, axes=(0,0)).transpose().astype('float32')


# get top sources for consideration
stdSources = (C.std(axis=0))*(A.mean(axis=(-1,-2)))
topRanks = np.argsort(stdSources)[-nSources:]


# create quantile-based threshold masks for top sources
topPixels = (A[topRanks].transpose(0,2,1)).reshape(A[topRanks].shape[0],-1)
maskQuantiles = np.array([ np.percentile(x[x>0], maskThresh) for x in topPixels] ).reshape((len(topRanks),1,1))
masks = A[topRanks].transpose(0,2,1) > maskQuantiles

# get df/f traces
df_f = cnm.estimates.F_dff

# get spikes
S = cnm.estimates.S

# get background
bg = np.tensordot(f.std(axis=0), b, axes=(0,0)).transpose().astype('float32')

saveDict = {'masks':masks,
            'df_f':df_f,
            'bg':bg,
            'A':A,
            'C':C,
            'b':b,
            'f':f,
            'S':S,
            'topRanks':topRanks
            }

savemat(outMatlabFile, saveDict)
np.savez(outArrayFile, **saveDict)


tif.imsave(os.path.join(outImageFolder, 'meanImageSources.tiff'), meanImageSources, imagej=True)
tif.imsave(os.path.join(outImageFolder, 'meanImage.tiff'), meanImage, imagej=True)
tif.imsave(os.path.join(outImageFolder, 'stdImageSources.tiff'), stdImageSources, imagej=True)
tif.imsave(os.path.join(outImageFolder, 'stdImage.tiff'), stdImage, imagej=True)

for i in range(len(masks)):
    img = Image.fromarray(255*masks[i].astype('uint8'))
    imgPath = os.path.join(outImageFolder, 'mask{0:0=03d}.png'.format(i))
    img.save(imgPath)

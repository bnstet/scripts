"""
Generate Df/f data from caiman .hdf5 outputs and motion-corrected movie files 

command format: genDff.py mcMovList.txt cmFileList.txt

files are saved in the cmFileList source folders as both npz and mat 
"""
import bokeh.plotting as bpl
import cv2
import glob
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import tifffile as tif
from time import time, sleep


try:
    cv2.setNumThreads(0)
except():
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
from caiman.utils.visualization import plot_contours, nb_view_patches, nb_plot_contour

from caiman.source_extraction import cnmf as cnmf
from scipy.sparse import csc_matrix
from scipy.misc import imresize
from scipy.signal import welch, fftconvolve
from scipy.interpolate import UnivariateSpline

from sklearn.decomposition import PCA
from sklearn import linear_model
from sklearn.metrics import r2_score

from scipy.ndimage.filters import gaussian_filter
from scipy.io import savemat

from caiman.source_extraction.cnmf.utilities import extract_DF_F, detrend_df_f

import pykalman as pk 
from matplotlib.colors import LogNorm

import h5py
bpl.output_notebook()

trialSize=60



###########

def fitTraceGLM(inTrace, templateTrace):
	inTrace = inTrace.reshape(-1,1)
	templateTrace = templateTrace.reshape(-1,1)
	reg = linear_model.LinearRegression()
	reg.fit(inTrace, templateTrace)
	return (reg.predict(inTrace).reshape((-1,)), reg)

def extractMaskAndBgTrace(mask, origMov, excludeMaskInBg=True, totalMask=None):
	mask = mask / np.linalg.norm(mask)
	if excludeMaskInBg:
		if totalMask is not None:
			bgMask = (totalMask == 0)
		else:
			bgMask = (mask == 0)
		bgMask = bgMask / np.linalg.norm(bgMask)
		movMean = (origMov * bgMask.reshape((1,) + mask.shape ) ).sum(axis=(-2,-1))
	else:
		movMean = origMov.mean(axis=(-2,-1))
	maskTrace = (origMov * mask.reshape((1,) + mask.shape ) ).sum(axis=(-2,-1))
	bgTrace, reg = fitTraceGLM(movMean,maskTrace)
	return (maskTrace, bgTrace, reg)


def cellInfoCaimanHdf5(hdf5File):
	"""
	Given a CaImAn output .hdf5 file, returns the spatial cell profiles, cell signal traces, and background profiles and traces
	Output:
	(cellProfs, signalTraces, bgProfs, bgTraces)
	cellProfs is a numCells x xDim x yDim array of spatial cell profiles
	signalTraces is a  tBins x numCells array of fluorescence traces from each cell
	bgProfs: nBGComps x xDim x yDim
	bgTraces: tBins x nBGComps
	"""
	cellProfs = []
	signalTraces = []

	with h5py.File(hdf5File, mode='r') as cFile:
		est = cFile['estimates']
		Ainfo = est['A']
		signalTraces = np.transpose( np.array(est['C']) )
		Adata = np.array(Ainfo['data'])
		Aindices = np.array(Ainfo['indices'])
		Aindptr = np.array(Ainfo['indptr'])
		Ashape = np.array(Ainfo['shape'])
		A = csc_matrix((Adata,Aindices,Aindptr) , shape=Ashape ).transpose()
		A = np.array(A.todense())
		cellProfs = A.reshape((A.shape[0],np.array(est['dims'])[0] ,-1))
		f = np.array(est['f']).transpose()
		b = np.array(est['b']).transpose()
		b = b.reshape(b.shape[0],cellProfs.shape[1],cellProfs.shape[2])
	return (cellProfs, signalTraces, b, f)



##########

with open('mc_files.txt', 'r') as f:
    mcFileList = [x.replace('\n','') for x in f.readlines()]

with open('cm_files.txt', 'r') as f:
    cmFileList = [x.replace('\n','') for x in f.readlines()]

	
	
for mcMovie, hdf5File in zip(mcFileList, cmFileList ):
	# run if loading batch cnmf output
	cnm = cnmf.cnmf.load_CNMF(hdf5File)
	if len(cnm.estimates.dims)==1:
		cnm.estimates.dims = cnm.estimates.dims*2

	# run if loading online cnmf output
	#cnm = cnmf.online_cnmf.load_OnlineCNMF(hdf5File)
	#if len(cnm.estimates.dims)==1:
	#    cnm.estimates.dims = cnm.estimates.dims*2

	# run to load in memmapped mc movie or mc tif directory
	if os.path.isdir(mcMovie):
		files = sorted(glob.glob(os.path.join(mcMovie, '*.tiff')) + glob.glob(os.path.join(mcMovie, '*.tif')))
		mcMov = tif.TiffSequence(files)
		mcMov = mcMov.asarray()
	elif os.path.splitext(mcMovie)=='tif':
		mcMov = tif.imread(mcMovie)
	elif os.path.splitext(mcMovie)=='mmap':
		Yr, T, dims  = cm.load_memmap(mcMovie)
		mcMov = Yr.reshape(dims + (T,)).transpose((2,0,1))
		
	if len(mcMov.shape)>=4:
		channelInd = np.argmin(mcMov.shape)
		mcMov = mcMov.mean(axis=channelInd)



	A, C, b, f = cellInfoCaimanHdf5(hdf5File)

	# get trial average of array.
	# time index should be in position 0.
	def trialAvg(arr, trialSize, idxStart=0, idxEnd=None):
		if idxEnd is None:
			idxEnd = len(arr)
		inputLength = idxEnd - idxStart
		if (inputLength % trialSize) != 0:
			raise Exception('Input length {} incompatible with trial size {} (remainder {})'.format(inputLength, trialSize, inputLength % trialSize))
		arr = arr[idxStart:idxEnd]
		return arr.reshape((int(arr.shape[0]/trialSize),trialSize) + arr.shape[1:], order='C').mean(axis=0)

	# generate masks from caiman regions
	maskThresh= 50 #remove pixels below this percentile from masks
	AReshape = A.transpose(0,2,1).reshape(A.shape[0],-1)
	rawMasks = A.transpose(0,2,1)
	maskQuantiles = np.array([ np.percentile(x[x>0], maskThresh) for x in AReshape] ).reshape((len(A),1,1))
	masks = A.transpose(0,2,1) > maskQuantiles
	normedMasks = masks / np.linalg.norm(masks,axis=(-1,-2), keepdims=True)


	totMask = (normedMasks.sum(axis=0) > 0).astype('float32')
	totMask = totMask / np.linalg.norm(totMask)




	maskResiduals = np.zeros( (len(rawMasks), len(mcMov) - (len(mcMov) %trialSize) ) )
	maskTraces = maskResiduals.copy()
	fitBgTraces = maskResiduals.copy()
	fitR2s = np.zeros(len(rawMasks))
	for n,mask in enumerate(rawMasks):
		maskTraces[n], fitBgTraces[n], _ = extractMaskAndBgTrace(rawMasks[n],mcMov[len(mcMov)%trialSize:], excludeMaskInBg=True, totalMask=totMask)
		maskResiduals[n] = maskTraces[n] - fitBgTraces[n]
		fitR2s[n] = r2_score(maskTraces[n], fitBgTraces[n])

	cmDff = cnm.estimates.F_dff[:,cnm.estimates.F_dff.shape[-1] % trialSize:]
	estDff = maskResiduals/maskTraces.mean(axis=-1,keepdims=True)


	# save results
	savePathNPZ = os.path.join( os.path.dirname(hdf5File), os.path.basename(hdf5File).replace('caiman', 'dffStats').replace('hdf5', 'npz'))
	savePathMat = os.path.join( os.path.dirname(hdf5File), os.path.basename(hdf5File).replace('caiman', 'dffStats').replace('hdf5', 'mat'))
	saveDict = {
		'rawMasks':rawMasks,
		'totMask':totMask,
		'normedMasks':normedMasks,
		'maskResiduals':maskResiduals,
		'maskTraces':maskTraces,
		'fitBgTraces':fitBgTraces,
		'fitR2s':fitR2s,
		'cmDff':cmDff,
		'estDff':estDff
		}
		

	np.savez_compressed(savePathNPZ, **saveDict)
	savemat(savePathMat, saveDict, do_compression=True)





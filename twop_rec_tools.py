#!/usr/bin/env python
"""Miscellaneous functions for dealing with 2p recordings"""

import numpy as np
import tifffile as tf
from scipy.io import loadmat, savemat, whosmat
import json
import os
from pathlib import Path
import h5py
from scipy.sparse import csc_matrix

def cellInfoS2pMat(matFile, excludeNonCell=True):
    """Take in a path to a .mat-file object from suite2p and extract the cell footprints and signal traces.
    Input:`
    matFile:.mat object filepath
    vidDims: (tBins, xDim, yDim) dimensions of movie
    excludeNonCell: exclude "cells" that have stat.iscell==False in suite2p
    Output:
    (cellMasks, signalTraces)
    cellMasks is a numCells x xDim x yDim array of cell binary masks
    signalTraces is a  tBins x numCells array of fluorescence traces from each cell"""
    

    matLoad = loadmat(matFile)
    matStat = matLoad['stat'][0]
    matOps = matLoad['ops'][0]

    statFields = list(matStat.dtype.fields.keys())
    opsFields = list(matOps.dtype.fields.keys())

    matOps = matOps[0]


    tBins, xDim, yDim = [ matOps[opsFields.index(x)][0,0] for x in ['Nframes','Lx', 'Ly']]


    # extract cell masks
    isCell = np.array([x[statFields.index('iscell')][0,0] for x in matStat], dtype='bool')
    if excludeNonCell:
        matStat = matStat[isCell]
    nCells = len(matStat)
    cellMasks = np.zeros((nCells, xDim, yDim), dtype='uint8')
    for nCell in range(nCells):
        xInds = matStat[nCell][statFields.index('xpix')]
        yInds = matStat[nCell][statFields.index('ypix')]
        cellMasks[nCell, xInds, yInds] = 1
    
    # extract cell fluorescence traces
    matFcell = np.transpose(matLoad['Fcell'][0,0])
    if excludeNonCell:
        matFcell = matFcell[:,isCell]
    return (cellMasks, matFcell)
    
def loadNeurofinderRegions(nfFolder):
    """
    Load in region masks from a Neurofinder folder path.
    Output:
    nCells x xDim x yDim binary masks 
    """
    folderPath = Path(nfFolder)
    jsonFile = str(list(folderPath.rglob('regions.json'))[0])
    imgDims = tf.imread( str(next(folderPath.joinpath('images').glob('*.tiff'))) ).shape
    with open(jsonFile,'r') as f:
        regions = json.load(f)

    masks = np.zeros( (len(regions),) + imgDims)
    for n,s in enumerate(regions):
        coords = s['coordinates'] 
        masks[n][tuple(zip(*coords))] = 1

    return masks


def cellInfoCaimanHdf5(hdf5File):
    """
    Given a CaImAn output .hdf5 file, returns the spatial cell profiles and cell signal traces.
    Output:
    (cellProfs, signalTraces)
    cellProfs is a numCells x xDim x yDim array of spatial cell profiles
    signalTraces is a  tBins x numCells array of fluorescence traces from each cell
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

    return (cellProfs, signalTraces)



def imageIou(img1,img2):
    """
    Return the intersection-over-union (IOU) score of two binary mask images (numpy arrays)
    The images don't need to be binary - any nonzero entries will be interpreted as 1s. For this to work however there should be no negative entries.
    Output: IOU score (float , between 0.0 and 1.0 inclusive)
    """

    intersection = np.sum( img1 * img2 != 0)
    union = np.sum( (img1 + img2) != 0 )
    return intersection / union


def imageCorr(img1,img2):
    """
    Return the (Pearson) correlation coefficient between two images (numpy arrays)
    Output: correlation coefficient (float, between -1.0 and 1.0 inclusive)
    """
    return  np.corrcoef( img1.flatten(), img2.flatten() )[0,1]
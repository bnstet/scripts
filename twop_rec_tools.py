#!/usr/bin/env python
"""Miscellaneous functions for dealing with 2p recordings"""

import numpy as np
import tifffile as tf
from scipy.io import loadmat, savemat, whosmat
import json


def cellInfoFromMat(matLoad, vidDims, excludeNonCell=True):
    """Take in an already-loaded .mat-file object from suite2p and extract the cell footprints and signal traces.
    Input:`
    matLoad:.mat object loaded into python from suite2p output
    vidDims: (tBins, xDim, yDim) dimensions of movie
    excludeNonCell: exclude "cells" that have stat.iscell==False in suite2p
    Output:
    (cellMasks, signalTraces)
    cellMasks is a numCells x xDim x yDim array of cell binary masks
    signalTraces is a numCells x tBins array of fluorescence traces from each cell"""
    
    tBins, xDim, yDim = vidDims
    
    # extract cell masks
    matStat = matLoad['stat'][0]
    isCell = np.array([x[27][0,0] for x in matStat], dtype='bool')
    if excludeNonCell:
        matStat = matStat[isCell]
    nCells = len(matStat)
    cellMasks = np.zeros((nCells, xDim, yDim), dtype='uint8')
    for nCell in range(nCells):
        xInds = matStat[nCell][2]
        yInds = matStat[nCell][1]
        cellMasks[nCell, xInds, yInds] = 1
    
    # extract cell fluorescence traces
    matFcell = matLoad['Fcell'][0,0]
    if excludeNonCell:
        matFcell = matFcell[isCell]
    return (cellMasks, matFcell)
    
def loadNeurofinderRegions(jsonFile, imgDim):
    """
    Load in region masks from a json Neurofinder regions.json file.
    imgDim: (xDim, yDim) dimensions of images
    Output:
    nCells x xDim x yDim binary masks 
    """
    with open(jsonFile,'r') as f:
        regions = json.load(f)

    masks = np.zeros( (len(regions),) + imgDim )
    for n,s in enumerate(regions):
        coords = s['coordinates'] 
        masks[n][list(zip(*coords))] = 1

    return masks

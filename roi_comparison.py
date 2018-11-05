#!/usr/bin/env python
"""
Run to compare ROI results for Caiman and Suite2p to ground truth
"""

import twop_rec_tools as tp
import os
from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat, whosmat

from pandas import DataFrame

"""
neurofinderList = ['/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.00/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.01/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.02/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.03/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.01.00/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.02.00/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.03.00/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.04.00/']


caimanList = ['/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.00/images/caiman_output.npz',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.01/images/caiman_output.npz',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.02/images/caiman_output.npz',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.03/images/caiman_output.npz',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.01.00/images/caiman_output.npz',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.02.00/images/caiman_output.npz',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.03.00/images/caiman_output.npz',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.04.00/images/caiman_output.npz']



suite2pList = ['/gpfs/data/shohamlab/ben/segmentation_project/Segmentation_Packages/suite2p_settings/results/237405/0/F_237405_0_plane1.mat',
'/gpfs/data/shohamlab/ben/segmentation_project/Segmentation_Packages/suite2p_settings/results/237405/1/F_237405_1_plane1.mat',
'/gpfs/data/shohamlab/ben/segmentation_project/Segmentation_Packages/suite2p_settings/results/237405/2/F_237405_2_plane1.mat',
'/gpfs/data/shohamlab/ben/segmentation_project/Segmentation_Packages/suite2p_settings/results/237405/3/F_237405_3_plane1.mat',
'/gpfs/data/shohamlab/ben/segmentation_project/Segmentation_Packages/suite2p_settings/results/237405/4/F_237405_4_plane1.mat',
'/gpfs/data/shohamlab/ben/segmentation_project/Segmentation_Packages/suite2p_settings/results/237405/5/F_237405_5_plane1.mat',
'/gpfs/data/shohamlab/ben/segmentation_project/Segmentation_Packages/suite2p_settings/results/237405/6/F_237405_6_plane1.mat']
"""

#% for local testing
neurofinderList = [r'C:\Users\bnste\Documents\segmentation_project\neurofinder\data\neurofinder.00.02']
caimanList = [r'C:\Users\bnste\Documents\segmentation_project\neurofinder\data\neurofinder.00.02\results\results.hdf5']
suite2pList = [r'C:\Users\bnste\Documents\segmentation_project\neurofinder\data\neurofinder.00.02\results\F_237405_2_plane1.mat']


caimanCells = [tp.cellInfoCaimanHdf5(x)[0] for x in caimanList]
suite2pCells = [tp.cellInfoS2pMat(x)[0] for x in suite2pList]
neurofinderCells = [tp.loadNeurofinderRegions(x) for x in neurofinderList]


#% construct mean images and compute IOU score between neurofinder (ground truth) and estimation methods.
#% caiman mean image is converted to binary mask by taking the 90th percentile of nonzero entries.

caimanThresh = .9
caimanMeans = [ np.mean(x, axis=0) for x in caimanCells  ]
caimanMeanMask = [ (x > np.percentile(x[x.nonzero()], caimanThresh )).astype('uint8') for x in caimanMeans  ]
suite2pMeans = [ np.mean(x, axis=0) for x in suite2pCells  ]
suite2pMeanMask = [ (x > 0).astype('uint8') for x in suite2pMeans   ]
neurofinderMeans = [ np.mean(x, axis=0) for x in neurofinderCells  ]
neurofinderMeanMask = [ (x > 0).astype('uint8') for x in neurofinderMeans  ]


scores = [ { 'caiman/S2p/IOU': tp.imageIou(caimanMeanMask[i], suite2pMeanMask[i] ),
                'caiman/NF/IOU': tp.imageIou(caimanMeanMask[i], neurofinderMeanMask[i] ),
                'S2p/NF/IOU': tp.imageIou(neurofinderMeanMask[i], suite2pMeanMask[i] ),
                'caiman/S2p/Corr': tp.imageCorr(caimanMeans[i], suite2pMeans[i] ),
                'caiman/NF/Corr': tp.imageCorr(caimanMeans[i], neurofinderMeans[i] ),
                'S2p/NF/Corr': tp.imageCorr(neurofinderMeans[i], suite2pMeans[i] )
                }  for i in range(len(neurofinderMeans))]

for n,x in enumerate(scores):
        print('Scoring info for dataset {}: '.format(n))
        print('Caiman data location: {}'.format(caimanList[n]))
        print('Suite2p data location: {}'.format(suite2pList[n]))
        print('Neurofinder data location: {}'.format(neurofinderList[n]))
        for key,val in x.items():
            print("Method: {}  Score: {:8.4f}".format(key, val))
        print('\n')

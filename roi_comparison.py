#!/usr/bin/env python
"""
Run to compare ROI results for Caiman and Suite2p to ground truth
"""

import twop_rec_tools as tp
import os
from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat, whosmat


neurofinderList = ['/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.00/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.01/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.02/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.00.03/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.01.00/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.02.00/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.03.00/',
'/gpfs/data/shohamlab/ben/segmentation_project/neurofinder/data/neurofinder.04.00/']
neurofinderList = [str(list(Path(x).rglob('regions.json'))[0]) for x in neurofinderList ]


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


#% for local testing
neurofinderList = [r'C:\Users\bnste\Documents\segmentation_project\neurofinder\data\neurofinder.00.02']
neurofinderList = [str(list(Path(x).rglob('regions.json'))[0]) for x in neurofinderList ]
caimanList = [r'C:\Users\bnste\Documents\segmentation_project\neurofinder\data\neurofinder.00.02\results\results.hdf5']
suite2pList = [r'C:\Users\bnste\Documents\segmentation_project\neurofinder\data\neurofinder.00.02\results\F_237405_2_plane1.mat']


caimanCells = [tp.cellInfoCaimanHdf5(x)[0] for x in caimanList]
suite2pCells = [tp.loadNeurofinderRegions(x)[0] for x in suite2pList]
neurofinderCells = [tp.loadNeurofinderRegions for x in neurofinderList]

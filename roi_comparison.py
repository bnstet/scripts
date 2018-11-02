#!/usr/bin/env python
"""
Run to compare ROI results for Caiman and Suite2p to ground truth
"""

import twop_rec_tools as tp
import os
from pathlib import Path
import numpy as np

vidDims = (3024, 520, 520)

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



suite2pList = [0,1,2,3,4,5,6,7]
s2pBasePath = '/gpfs/data/shohamlab/ben/segmentation_project/Segmentation_Packages/suite2p_settings/results/237405'
suite2pList = [ str(list(Path(os.path.join(s2pBasePath,str(x))).rglob('*plane1.npz'))[0]) for x in suite2pList]


caimanDataList = [np.load(x) for x in caimanList]
suite2pDataList = [np.load(x) for x in suite2pList]


caimanCells = [x['A'] for x in caimanDataList]
suite2pCells = [tp.cellInfoFromMat(x, vidDims)[0] for x in suite2pDataList]
neurofinderCells = [tp.loadNeurofinderRegions(x, vidDims[1:]) for x in neurofinderList]

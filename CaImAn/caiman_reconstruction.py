#!/usr/bin/env python
"""
Run to examine CaImAn output (cell footprints, cell output traces, and more).
Reconstruct denoised movie from CaImAn output.
"""

import twop_rec_tools as tp
import os
from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat, whosmat

from pandas import DataFrame




caimanList = ['/gpfs/home/stetlb01/scripts/CaImAn/results/caiman_analysis_0.hdf5',
'/gpfs/home/stetlb01/scripts/CaImAn/results/caiman_analysis_1.hdf5',
'/gpfs/home/stetlb01/scripts/CaImAn/results/caiman_analysis_2.hdf5',
'/gpfs/home/stetlb01/scripts/CaImAn/results/caiman_analysis_3.hdf5',
'/gpfs/home/stetlb01/scripts/CaImAn/results/caiman_analysis_4.hdf5',
'/gpfs/home/stetlb01/scripts/CaImAn/results/caiman_analysis_5.hdf5',
'/gpfs/home/stetlb01/scripts/CaImAn/results/caiman_analysis_6.hdf5',
'/gpfs/home/stetlb01/scripts/CaImAn/results/caiman_analysis_7.hdf5',]



caimanInfo = [tp.cellInfoCaimanHdf5(x) for x in caimanList]
caimanA = [x[0] for x in caimanInfo]
caimanC = [x[1] for x in caimanInfo]


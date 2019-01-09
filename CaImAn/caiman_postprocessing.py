import numpy as numpy
import caiman as cm
import tifffile as tif

import os
import sys
import glob

import twop_rec_tools as tp

import argparse



parser = argparse.ArgumentParser(description='Utilize CaImAn hdf5 analysis output to produce summary images and masks')
parser.add_argument('infile', type=str, nargs=1, help='input hdf5 file')
parser.add_argument('outdir',type=str, nargs=1, help='directory in which to store outputs')

args = parser.parse_args()

infile = args.infile[0]
outdir = args.outdir[0]

print(infile, outdir)
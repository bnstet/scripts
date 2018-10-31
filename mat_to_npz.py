#!/usr/bin/env python

"""
Script for converting Matlab .mat files into numpy .npz files. 

Usage:
./mat_to_npz.py infile [-o outfile] [-c]
./mat_to_npz.py dir [-r] [-c]

If called on a file, will only convert that file. If the output name is not explicitly given with -o, then the output file name will be the original name with .mat replaced by .npz.

If called on a directory, all .mat files within will be converted to .npz files with default names. The -r flag searches recursively within the directory; all output files will be located on the level of the corresponding input file.

The -c flag causes the output to be saved in compressed .npz format; the default save method is noncompressed.
"""

from __future__ import print_function

###

from numpy import savez, savez_compressed
from scipy.io import loadmat, whosmat
from pathlib import Path
import argparse
import os
import sys

###

parser = argparse.ArgumentParser()
parser.add_argument('fileOrDir', type=str, help='Input .mat file, or directory in which to look for .mat files')
parser.add_argument('-o', type=str, help='output file name (optional)', default=None)
parser.add_argument('-r', action='store_true', help='recursively search for .mat files in directory')
parser.add_argument('-c', action='store_true', help='save .npz file(s) in compressed format')
args = parser.parse_args();

###

def mat_to_npz(inFile, outFile=None, comp=False):
    mat_contents = loadmat(inFile)
    saveFn = savez_compressed if comp else savez
    outFile = outFile if not (outFile is None) else inFile.replace('.mat','.npz')
    saveFn(outFile, **mat_contents)



###

isDir = os.path.isdir(args.fileOrDir)
isFile = os.path.isfile(args.fileOrDir)




if (not isDir) and (not isFile):
    sys.exit('Provided file or directory does not exist.')
if isDir:
    if args.r:
        fileList = list(Path(args.fileOrDir).rglob("*.mat"))
    else:
        fileList = list(Path(args.fileOrDir).glob("*.mat"))
    print('Converting {} .mat files in {}'.format(len(fileList), args.fileOrDir))
    for n,myFile in enumerate(fileList,1):
        print('({}/{}) Converting file {}'.format(n,len(fileList),myFile))
        mat_to_npz(str(myFile), outFile=None, comp=args.c)
else:
    mat_to_npz(fileOrDir, outFile=args.outFile, comp = args.c)

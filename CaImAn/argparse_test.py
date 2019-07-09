import os,sys,argparse

parser = argparse.ArgumentParser(description='Execute the CaImAn (Calcium Imaging Analysis) pipeline on a video recording (.tif file or .tiff stack). Motion correct, determine neuronal regions, extract signals and denoise.')

parser.add_argument('--strvar', type=str, nargs=1, help='file name (for .tif file) or folder (for .tiff stack). can also be a .txt file containing a list of .tif files for multiple chunks.')
parser.add_argument('--strlistvar',type=str, nargs='?', default=None, help='template .tif file for motion correction')
parser.add_argument('--strlistvar2',type=str, nargs='2', default=None, help='template .tif file for motion correction')
parser.add_argument('--boolvar', action='store_true', help='if used, then no motion correction will be run on the input video/mmap file')
arser.add_argument('--intvar', type=int, nargs='1', default=0, help='slurm ID for multi-node processing. leave out to default to 0')
parser.add_argument('--intlistvar', type=int, nargs='?', default=0, help='slurm ID for multi-node processing. leave out to default to 0')
parser.add_argument('--intlistvar2', type=int, nargs='2', default=(1,2), help='slurm ID for multi-node processing. leave out to default to 0')


args = vars(parser.parse_args())

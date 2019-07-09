
import os
import caiman as cm
import argparse

parser = argparse.ArgumentParser(description='save multiple video files as a single tif')
parser.add_argument('infile', type=str, nargs=1, help='.txt file containing chunk file locations')
parser.add_argument('outfile', type=str, nargs=1, help='output file name (.tif file name should be used)')

args = vars(parser.parse_args())
infile = args['infile'][0]
outfile = args['outfile'][0]


if '.txt' in infile:
    with open(infile, 'r') as f:
        flist = []
        for line in f:
            flist = flist + [line.replace('\n','')]


m = cm.load_movie_chain(flist)
m.save(outfile)

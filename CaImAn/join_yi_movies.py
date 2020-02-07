from caiman import movie, load_movie_chain
import tifffile as tif, numpy as np, os, glob

import argparse, warnings


parser = argparse.ArgumentParser(description='join together a directory of movies with yi\'s tagging system, in order, into a combined tif file')

parser.add_argument('indir', type=str, nargs=1, help='directory name containing movie files or movie stack directories')
parser.add_argument('--match', type=str, nargs='?', default=None, help='matching pattern for movie names (default is "_mc" or "_mc.tif")')
parser.add_argument('--movlength', type=int, nargs='?', default=None, help='use to cut off the first frames of constituent movies to force this length')



args=parser.parse_args()

if isinstance(args.indir, type([])):
    args.indir = args.indir[0]


if args.match is None:
    mov_files = glob.glob(os.path.join(args.indir, "*_mc")) + glob.glob(os.path.join(args.indir, "*_mc.tif"))
else:
    mov_files = glob.glob(os.path.join(args.indir, args.match))






if len(mov_files)==0:
    raise FileNotFoundError('Cannot find matching movie files in directory {}'.format(args.indir))


mov_files = sorted(mov_files, key = lambda x : int(x.split('-')[-1].split('_')[0]))
print('sorted movie list:')
[print(x) for x in mov_files]



def mov_load_func(x):
    with warnings.catch_warnings():
        if os.path.isdir(x):
            return load_movie_chain(sorted(glob.glob(os.path.join(x,'*.tif*'))))
        else:
            return load_movie_chain(x)


mov_list = [mov_load_func(x) for x in mov_files]


if args.movlength is not None:
    for n,mov in enumerate(mov_list):
        if mov.shape[0] < args.movlength:
            raise Exception('Movie {} has length {} < {}'.format(mov_files[n], mov.shape[0], args.movlength))
    mov_list = [x[-args.movlength:] for x in mov_list]

mov = movie(np.concatenate(mov_list, axis=0))


saveloc = os.path.join(args.indir,'combined.tif')
mov.save(saveloc)
print('combined movie saved to {}'.format(saveloc)) 

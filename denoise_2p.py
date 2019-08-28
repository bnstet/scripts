"""
Script to process an input 2p video from a .tif(f) file and denoise it.

Default output is a .tif file containing the denoised movie.

python denoise_2p.py file_name1 [filename2] ... [--save_npz] [--npz_only] [--out_res]


--save_npz [file_name]: save the outputs of the denoising/compression algorithm to file_name as an npz file
--npz_only: if used, do not output denoised tif, just the npz file. can only be used with save_npz.
--out_res: change output resolution to this pair. this is applied before any denoising.
--block_ratio: scale down each side for the algorithm block decompositon. example input: --block_ratio 4 5 would result in a 20 block decomposition.
--max_comp: maximum components per block
--no_rescale: if flagged, do not rescale the output to return it to the original intensity distribution; instead leave it in the "spatially uniform noise" state
--pca_comp: number of components to use for the preliminary pca decomposition
"""

# handle input arguments
parser = argparse.ArgumentParser(description='filter, smooth, and resize 2p video recordings')
parser.add_argument('file_list', type=str, nargs='+', help='list of files. if one .txt file is provided, that file will be read to retrieve the contained file list (line-separated)')
parser.add_argument('--out_res', type=int, nargs='+', default=None, help='output resolution (defaults to preserving resolution)')
parser.add_argument('--n_channels', type=int, default=1, help='number of input channels')
parser.add_argument('--saved_channels', type=int, nargs='+', default=[0], help='input channels to save (0-indexed)')
parser.add_argument('--pct_thresh', type=float, default=99.99, help='pixel percentile threshold; use 100 for no percentile filtering')
parser.add_argument('--alpha', type=float, default=0, help='AR filter smoothing parameter, alpha=1 means a flat moving average of size n_smooth, while alpha=0 means no smoothing. default 0')
parser.add_argument('--n_smooth', type=int, default=1, help='size of smoothing window. 1 frame is the default')
parser.add_argument('--keep_bad_frames', action='store_true', help='if this flag is activated, keep outlier frames instead of interpolating through them')
parser.add_argument('--ds_factor', type=int, default=1, help='downsample factor; one frame is kept every ds_factor frames. applied prior to smoothing!')
parser.add_argument('--med_ds', dest='ds_func', action='store_const', const=np.median, default=np.mean, help='use median instead of mean for downsampling')


args = parser.parse_args()

print('Using arguments: {}'.format(args),flush=True)

out_res = tuple(args.out_res) if args.out_res is not None else None
n_channels = args.n_channels
saved_channels = args.saved_channels
filter_frames = not args.keep_bad_frames
pct_thresh = args.pct_thresh
file_list = args.file_list
alpha = args.alpha
nSmooth = args.n_smooth
ds_factor = args.ds_factor
ds_func = args.ds_func

# get file list from text file if text file is in the file_list input variable
if '.txt' in file_list[0]:
    with open(file_list[0], mode='r') as f:
        file_list = [x.replace('\n','') for x in f.readlines()]

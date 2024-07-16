import argparse
from pathlib import Path
from astropy.stats import mad_std
import ccdproc as ccdp
import numpy as np

#bias_dir = Path('bias_1_1')
#
#files = ccdp.ImageFileCollection(bias_dir)
#
#bias_files = files.files_filtered(imagetyp='BIAS', include_path=True)
#
#combined_bias = ccdp.combine(bias_files,method='average',
#                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
#                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
#
#combined_bias.meta['combined'] = True
#combined_bias.write(bias_dir / 'combined_bias_1_1.fit', overwrite=True)
#
#
#bias_dir = Path('bias_2_2')
#
#files = ccdp.ImageFileCollection(bias_dir)
#
#bias_files = files.files_filtered(imagetyp='BIAS', include_path=True)
#
#combined_bias = ccdp.combine(bias_files,method='average',
#                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
#                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
#
#combined_bias.meta['combined'] = True
#combined_bias.write(bias_dir / 'combined_bias_2_2.fit', overwrite=True)


def process_bias_directory(bias_dir):
    files = ccdp.ImageFileCollection(bias_dir)
    bias_files = files.files_filtered(imagetyp='BIAS', include_path=True)
    
    combined_bias = ccdp.combine(bias_files, method='average',
                                 sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                 sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
    
    combined_bias.meta['combined'] = True
    combined_bias.write(bias_dir / 'combined_bias.fit', overwrite=True)

def main():
    parser = argparse.ArgumentParser(description="Process bias files in the given directories.")
    parser.add_argument('directories', metavar='DIR', type=str, nargs='+',
                        help='one or more directories containing bias files')

    args = parser.parse_args()

    for directory in args.directories:
        bias_dir = Path(directory)
        if bias_dir.is_dir():
            process_bias_directory(bias_dir)
        else:
            print(f"Error: {bias_dir} is not a valid directory")

if __name__ == "__main__":
    main()

import os
import argparse
#from astropy.nddata import CCDData
from ccdproc import CCDData, Combiner, subtract_bias
import numpy as np
from astropy.io import fits
from astropy import units as u

def inv_median(a):
    return 1 / np.median(a)

def get_header_values(fits_files, keyword):
    """
    Collects values of a specified header keyword from a list of FITS files.

    Parameters:
        fits_files (list): List of paths to FITS files.
        keyword (str): Header keyword to retrieve.

    Returns:
        list: Values of the header keyword from each FITS file.
    """
    header_values = []
    for file_path in fits_files:
        with fits.open(file_path) as hdul:
            header = hdul[0].header
            value = header.get(keyword, 'N/A')  # Default to 'N/A' if keyword does not exist
            header_values.append(value)
    return header_values

def create_flat_field_image(file_paths, header_values, bias, outdir, unit='adu'):
    """
    Creates a flat field image by averaging the data from provided FITS files using ccdproc.

    Parameters:
        file_paths (list): List of FITS file paths to combine.

    Returns:
        CCDData: Combined flat field image.
    """
    index_dict = {}
    for index, angle in enumerate(header_values):
        if angle not in index_dict:
            index_dict[angle] = []
        index_dict[angle].append(index)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for value, indices in index_dict.items():
        ccds = [CCDData.read(file_paths[i], unit=unit) for i in indices]
#        ccds = [subtract_bias(ccd, bias) for ccd in ccds]

        # Create a circular mask
        masks = []
        for ccd in ccds:
            y, x = np.indices(ccd.shape)
            center_x, center_y = ccd.shape[1] / 2, ccd.shape[0] / 2
            mask = ((x - 200)**2 + (y - 206)**2) <= 130**2
            masks.append(mask)

        combiner = Combiner(ccds)

        # Apply the mask and calculate scaling using only the masked region
        combiner.scaling = [1./np.mean(ccd.data[mask]) for ccd, mask in zip(ccds, masks)]
        # combiner.scaling = [1./np.median(ccd.data) for ccd in ccds]

        # clip
        combiner.sigma_clipping(low_thresh=3, high_thresh=3, func='median', dev_func='std')
        
        combined_image = combiner.average_combine()

        #combined_image = combiner.average_combine(scale=inv_median)
        file_name = os.path.join(outdir, value.replace(' ', '_') + '.fits')
        combined_image.write(file_name, overwrite=True)

def main():
    parser = argparse.ArgumentParser(description='Extract header information and create a flat field image from FITS files using ccdproc.')
    parser.add_argument('fits_files', nargs='+', help='Paths to FITS files')
    parser.add_argument('--keyword', type=str, default='OBJECT', help='Header keyword to extract (default: OBJECT)')
    parser.add_argument('--bias_path', required=True, help="Path to the master bias FITS file.")
    parser.add_argument('--outdir', type=str, help='Output FITS file path')

    args = parser.parse_args()

    # read in master bias
    master_bias = CCDData.read(args.bias_path, unit='adu')
    
    header_values = get_header_values(args.fits_files, args.keyword)
    #print("Header values extracted:", header_values)

    flat_field_image = create_flat_field_image(args.fits_files,header_values,master_bias,args.outdir)

if __name__ == '__main__':
    main()

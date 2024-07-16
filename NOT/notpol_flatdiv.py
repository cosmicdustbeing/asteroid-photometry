import os
import argparse
from astropy.io import fits
from ccdproc import CCDData, flat_correct, subtract_bias
import numpy as np

def find_flat_filename(object_header, flat_files):
    """
    Extracts the angle from the object header and matches it to the corresponding flat filename.
    """
    angle = object_header.split()[-1]  # Assumes the angle is the last part of the header string
    for flat_file in flat_files:
        if f"_{angle}.fits" in flat_file:
            return flat_file
    return None

def bin_image(image, binning_factor):
    """
    Bins the image by the given binning factor.
    
    Parameters:
        image (ndarray): The input image to be binned.
        binning_factor (int): The factor by which to bin the image (e.g., 2 for 2x2 binning).
    
    Returns:
        ndarray: The binned image.
    """
    shape = (image.shape[0] // binning_factor, binning_factor,
             image.shape[1] // binning_factor, binning_factor)
    binned_image = image.reshape(shape).mean(axis=(1, 3))
    return binned_image

def divide_by_flat_field(image_paths, bias_path, flat_dir, output_dir, unit='adu'):
    """
    Divides each image by its corresponding flat field image based on the angle mentioned in the OBJECT header.
    
    Parameters:
        image_paths (list): List of paths to the target FITS files.
        bias_path (str): Path to the master bias FITS file.
        flat_dir (str): Directory containing the flat field FITS files.
        output_dir (str): Directory to save the corrected FITS files.
        unit (str): Unit of the data in the FITS files (default is 'adu').
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # load the bias frame
    master_bias = CCDData.read(bias_path, unit=unit)
        
    # Load all flat files from the directory
    flat_files = [f for f in os.listdir(flat_dir) if f.endswith('.fits')]
    flat_files = [os.path.join(flat_dir, f) for f in flat_files]

    for image_path in image_paths:
        with fits.open(image_path) as hdul:
            object_header = hdul[0].header['OBJECT']
            detxbin = hdul[0].header.get('DETXBIN', 1)
            detybin = hdul[0].header.get('DETYBIN', 1)
            flat_filename = find_flat_filename(object_header, flat_files)
            
            if flat_filename:
                with fits.open(flat_filename) as flat_hdul:
                    # Convert both science image and flat field image to CCDData objects
                    ccd_image = CCDData(hdul[1].data, unit=unit, meta=hdul[0].header+hdul[1].header)

                    bias_subtracted = subtract_bias(ccd_image, master_bias)

                    if detxbin == 1 and detybin == 1:
                        binned_data = bin_image(bias_subtracted.data, 2)
                        binned_image = CCDData(binned_data, unit=unit, meta=bias_subtracted.meta)
                        binned_image.header['DETXBIN'] = 2
                        binned_image.header['DETYBIN'] = 2
                    else:
                        binned_image = bias_subtracted
                        
                    ccd_flat = CCDData(flat_hdul[0].data, unit=unit, meta=flat_hdul[0].header)
                    
                    # Perform flat field correction
                    corrected_ccd = flat_correct(binned_image, ccd_flat, min_value=0.1)

                    # Create a new FITS file with the corrected data
                    output_path = os.path.join(output_dir, os.path.splitext(os.path.basename(image_path))[0] + '_cal.fits')
                    corrected_ccd.write(output_path, overwrite=True)
                    print(f"Processed {image_path} with {flat_filename}, saved to {output_path}")
            else:
                print(f"No corresponding flat file found for {image_path}")

def main():
    parser = argparse.ArgumentParser(description="Divide FITS images by corresponding flat field images.")
    parser.add_argument('image_paths', nargs='+', help="Paths to the target FITS files.")
    parser.add_argument('--bias_path', required=True, help="Path to the master bias FITS file.")
    parser.add_argument('--flat_dir', required=True, help="Directory containing the flat field FITS files.")
    parser.add_argument('--output_dir', required=True, help="Directory to save the corrected FITS files.")
    parser.add_argument('--unit', default='adu', help="Unit of the data in the FITS files (default is 'adu').")
    
    args = parser.parse_args()
    
    divide_by_flat_field(args.image_paths, args.bias_path, args.flat_dir, args.output_dir, args.unit)

if __name__ == "__main__":
    main()

import os
import math
from astropy.io import fits
import warnings
from glob import glob
import photutils
from reproject import reproject_interp  # https://reproject.readthedocs.io/en/stable/
from tqdm import tqdm

class PrepImg:
    
    ### INITIALIZE
    ### Only inputs are image directory, file naming pattern, and field name

    def __init__(self, image_dir, file, field):
        # Initialize image locations and future names
        self.image_dir = image_dir
        self.file = file
        self.field = field
        self.imfile = f'{self.field}_*_sci.fits'

        # Create empty data structures
        self.i2d_img = []
        self.sci_img = []
        self.wht_img = []
        self.bkg_img = []
        
        # Ignore astropy FITSFixedWarning
        warnings.simplefilter('ignore')

        ### GET I2D IMAGES
        ### Uses image_dir and file to find existing i2d images

        self.i2d_img = glob(os.path.join(self.image_dir, self.file))
        self.i2d_img.sort()

        ### SPLIT I2D FILES
        ### Creates 'sci' and 'wht' files from the i2d images
        
        for input_image in tqdm(self.i2d_img,desc='Saving individual images...'):
            for i in range(len(input_image)):
                if input_image[i] == 'f' and (input_image[i+4] == 'w' or input_image[i+4] == 'm') and input_image[i+1] in '0 1 2 3 4 5 6 7 8 9'.split():
                    filt = input_image[i:i+5]
                    break
            
            hdu_list = fits.open(input_image)
            header = hdu_list['sci'].header[:]
            
            # Only the science extension header has any info including WCS
            # So we'll just use that for all the output files

            for extension in 'sci wht'.split():
                output_image = os.path.join(self.image_dir, f'{self.field}_{filt}_{extension}.fits')
                
                if extension == 'sci':
                    self.sci_img.append(output_image)
                else:
                    self.wht_img.append(output_image)
                
                header['EXTNAME'] = extension
                
                fits.writeto(output_image, hdu_list[extension].data, header, overwrite=True)
                
        ### GET ZEROPOINTS
        ### Use pixel size in image header to calculate photometric zeropoints for SW and LW filters
        ### If SW and LW images have the same pixel size, zp_sw == zp_lw

        zeropoints = []
        for file in tqdm(self.sci_img, desc="Calculating zeropoints..."):
            hdul = fits.open(file)
            header = hdul[0].header
                        
            pixel_size = header['PIXAR_A2']**0.5
            zp = 8.9 - 2.5 * math.log10(1e+6 / ( (360 * 3600) / (2 * math.pi * pixel_size) )**2)
            zeropoints.append(zp)
                
            hdul.close()
            
        self.zp_sw = max(zeropoints)
        self.zp_lw = min(zeropoints)
                
        ### REPROJECT SCI AND WHT FILES
        ### Having the images on the same pixel grid will be beneficial for SExtractor

        images = self.sci_img + self.wht_img
        for image in tqdm(images, desc='Reprojecting pixel grids...'):
            ref_image = self.sci_img[0]

            hdul = fits.open(ref_image)
            ref_header = hdul['sci'].header

            hdu = fits.open(image)
            data = hdu[0]

            reprojected_data, footprint = reproject_interp(data, ref_header)
            
            fits.writeto(image, reprojected_data.astype('float32'), ref_header, overwrite=True)
            hdu.close()

        ### DONE SPLITTING AND REPROJECTING
        ### Print
        
        print(f'Field: {self.field}\nImage Directory: {self.image_dir}\nImage File Pattern: {self.imfile}\n' +
              f'SW ZP: {self.zp_sw:.3f}\nLW ZP: {self.zp_lw:.3f}')
        
    ### Method bkgsub
    ### optionally subtract the background from the science images
    ### size is the size of each square that background is measured and subtracted from

    def bkgsub(self, size=100):
        for image in tqdm(self.sci_img, desc='Background Subtraction...'):
            hdul = fits.open(image)
            data = hdul['sci'].data

            background_map = photutils.Background2D(data, size, filter_size=5)  
            data = data - background_map.background.astype('float32')

            hdul['sci'].data = data
            newfile = image.replace(".fits", "_bkgsub.fits")
            self.bkg_img.append(newfile)
            hdul.writeto(newfile, overwrite=True)
            hdul.close()
            
        print(f'Field: {self.field}\nImage Directory: {self.image_dir}\nImage File Pattern: ' +
              f'{self.imfile.replace(".fits", "_bkgsub.fits")}')
import os
import math
from astropy.io import fits
from glob import glob
import photutils
from reproject import reproject_interp  # https://reproject.readthedocs.io/en/stable/

class prepimg:
    
    def __init__(self):
        self.image_dir = 'N/A'
        self.file = 'N/A'
        self.imfile = 'N/A'
        self.field = 'N/A'
        self.zp_sw = 28.0865
        self.zp_lw = 28.0865
        self.i2d_img = []
        self.sci_img = []
        self.wht_img = []
        self.bkg_img = []
        
    def get_i2d(self):
        self.i2d_img = glob(os.path.join(self.image_dir, self.file))
        
    def split_i2d(self):
        for input_image in self.i2d_img:
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
                
                print('SAVING', output_image)
                fits.writeto(output_image, hdu_list[extension].data, header, overwrite=True)
                
    def get_zp(self):
        zeropoints = []
        for file in self.sci_img:
            hdul = fits.open(file)
            header = hdul[0].header
                        
            pixel_size = header['PIXAR_A2']**0.5
            zp = 8.9 - 2.5 * math.log10(1e+6 / ( (360 * 3600) / (2 * math.pi * pixel_size) )**2)
            zeropoints.append(zp)
                
            hdul.close()
            
        self.zp_sw = max(zeropoints)
        self.zp_lw = min(zeropoints)
                
    def reproject_sci(self):
        for image in self.sci_img + self.wht_img:
            ref_image = self.sci_img[0]
            hdul = fits.open(ref_image)
            ref_header = hdul['sci'].header
            print("Reprojecting...")  # 1 minute
            hdu = fits.open(image)
            data = hdu[0]
            reprojected_data, footprint = reproject_interp(data, ref_header)
            fits.writeto(image, reprojected_data.astype('float32'), ref_header, overwrite=True)
            hdu.close()

    def go(self, image_dir, file, field):
        self.image_dir = image_dir
        self.file = file
        self.field = field
        self.imfile = f'{field}_*_sci.fits'
        
        self.get_i2d()
        self.split_i2d()
        self.get_zp()
        self.reproject_sci()
        
        print(f'Field: {self.field}\nImage Directory: {self.image_dir}\nImage File Pattern: {self.imfile}\n' +
              f'SW ZP: {self.zp_sw}\nLW ZP: {self.zp_lw}')
        
    def bkgsub(self, size=100):
        for image in self.sci_img:
            hdul = fits.open(image)
            data = hdul['sci'].data

            background_map = photutils.Background2D(data, size, filter_size=5)    
            data = data - background_map.background.astype('float32')

            hdul['sci'].data = data
            newfile = image.replace(".fits", "_bkgsub.fits")
            self.bkg_img.append(newfile)
            print('SAVING', newfile)
            hdul.writeto(newfile, overwrite=True)
            hdul.close()
            
        print(f'Field: {self.field}\nImage Directory: {self.image_dir}\nImage File Pattern: ' +
              f'{self.imfile.replace(".fits", "_bkgsub.fits")}')
        

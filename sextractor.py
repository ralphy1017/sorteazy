import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.table import Table, vstack
import math
import statistics
from glob import glob

class sextractor:
    def __init__(self): #make empty data structures and set default psf and sb dictionaries
        self.imfiles = {}
        self.whtfiles = {}
        self.areas = []
        self.exposures = []
        self.catnames = []
        self.initialize()
        self.auto_psf()
        self.auto_sb()
    
    def initialize(self, field='N/A', zp=28.0865, config_file='N/A', cat_dir='N/A', image_dir='N/A', imfile='N/A'):
        self.field = field #get field, zeropoint, SExtractor config file, catalog directory, and image directory
        self.zeropoint = zp
        self.config_file = config_file
        self.cat_dir = cat_dir
        self.image_dir = image_dir
        self.imfile = imfile
        if self.field != 'N/A':
            print(f'SExtractor initialized:\nField: {self.field}\nZeropoint: {self.zeropoint}\nConfig File:' +
                  f'{self.config_file}\nCatalog Directory: {self.cat_dir}\nImage Directory: {self.image_dir}' +
                  f'\nImage File Pattern: {self.imfile}')

    def auto_psf(self): #default psf dictionaries
        self.psf_fwhm = {'F090W': 0.06, 'F115W': 0.05, 'F150W': 0.055, 'F182M': 0.06, 'F200W': 0.065, 
                         'F210M': 0.07, 'F277W': 0.1, 'F300M': 0.11, 'F335M': 0.12, 'F356W': 0.15, 
                         'F360M': 0.15, 'F410M': 0.15, 'F444W': 0.15}
        self.psf_max = {'F090W': 0.11, 'F115W': 0.11, 'F150W': 0.08, 'F182M': 0.085, 'F200W': 0.09,
                        'F210M': 0.1, 'F277W': 0.13, 'F300M': 0.14, 'F335M': 0.15, 'F356W': 0.18,
                        'F360M': 0.17, 'F410M': 0.17, 'F444W': 0.17}
    
    def auto_sb(self): #default sensitivity dictionaries
        self.sb_min = {'F090W': 30, 'F115W': 29.3, 'F150W': 29.5, 'F182M': 28.1, 'F200W': 29.7,
                       'F210M': 27.8, 'F277W': 30.2, 'F300M': 29, 'F335M': 29, 'F356W': 30.1,
                       'F360M': 28.9, 'F410M': 29, 'F444W': 30}
        self.sb_turn = {'F090W': 0.2, 'F115W': 0.16, 'F150W': 0.21, 'F182M': 0.16, 'F200W': 0.19,
                        'F210M': 0.17, 'F277W': 0.2, 'F300M': 0.18, 'F335M': 0.18, 'F356W': 0.22,
                        'F360M': 0.18, 'F410M': 0.2, 'F444W': 0.2}
        self.sb_slope = {'F090W': -5, 'F115W': -4.8, 'F150W': -5.2, 'F182M': -5, 'F200W': -5.1,
                         'F210M': -5, 'F277W': -5, 'F300M': -5, 'F335M': -5, 'F356W': -4.9,
                         'F360M': -5, 'F410M': -4.9, 'F444W': -5}
            
    def set_psf(self, filt, fwhm, end): #minimum and maximum FWHM to classify an object as a star
        self.psf_fwhm[filt] = fwhm
        self.psf_max[filt] = end
        
    def set_sb(self, filt, minimum, turn, slope): #sets surface brightness / point source sensitivity curve
        self.sb_min[filt] = minimum
        self.sb_turn[filt] = turn
        self.sb_slope[filt] = slope
        
    def get_images(self): #get an invdividual image and weight file
        #imfile should be a filename pattern (e.g., 'jwst_field_*_30mas_sci.fits')
        #the filter should be replaced by * to get all files of that pattern
        #the file pattern should have a 'wht' file for each 'sci' file
        files = glob(os.path.join(self.image_dir, self.imfile))
        for image in files:
            imfile = image.lower().split('/')[-1]
            for i in range(len(imfile)):
                if imfile[i] == 'f' and (imfile[i+4] == 'w' or imfile[i+4] == 'm') and imfile[i+1] in '0 1 2 3 4 5 6 7 8 9'.split():
                    filt = imfile[i:i+5]
                    break
            self.imfiles[filt.upper()] = image
            self.whtfiles[filt.upper()] = image.replace('sci', 'wht').replace('SCI', 'WHT').replace('drz', 'wht')
        print('Images Found')
    
    def get_areas(self): #used for number counts
        areas = []
        for file in self.whtfiles.values():
            hdul = fits.open(file)
            data = hdul[0].data
            idx = np.where(data > 0)
            hdul.close()
            areas.append(len(idx[0]) * 0.0009 / 60**4)
        self.areas = areas
        print('Areas calculated.')
        
    def get_exposures(self): #input to source extractor as gain = exposure time
        for i in range(len(self.imfiles)):
            hdul = fits.open(list(self.imfiles.values())[i])
            header = hdul[0].header
            filt = list(self.imfiles.keys())[i]
            if int(filt[1:4]) < 250: #related to central wavelength (SW)
                self.exposures.append(header['XPOSURE'] / 8)
            elif int(filt[1:4]) >= 250: #related to central wavelength (LW)
                self.exposures.append(header['XPOSURE'] / 2)
            else:
                print('header matching failed')
            hdul.close()
        print('Exposure times found.')
    
    def sextract(self, overwrite=True): #run source extractor
        for i in range(len(self.imfiles)):
            catname = os.path.join(self.cat_dir, f'{self.field}_{list(self.imfiles.keys())[i]}_cat.txt')
            self.catnames.append(catname)
            if not os.path.exists(catname) or overwrite == True:
                os.system(f'sex {list(self.imfiles.values())[i]} -c {self.config_file} -MAG_ZEROPOINT {self.zeropoint} -CATALOG_NAME {catname} -GAIN {self.exposures[i]} -WEIGHT_IMAGE {list(self.whtfiles.values())[i]} -WEIGHT_TYPE MAP_WEIGHT')
        print('Source Extractor done.')
        
    def go(self): #run code after running self.initialize() and self.get_image()
        if self.field == 'N/A':
            print('Cannot go without initializing field (call sextractor.initialize())')
            
        self.get_images()
        self.get_areas()
        self.get_exposures()
        self.sextract()
        

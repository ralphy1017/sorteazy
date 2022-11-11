import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.table import Table
import math
from glob import glob
                    
class sextractor:
    
    def __init__(self): #make empty data structures and set default psf and sb dictionaries
        self.imfiles = {}
        self.whtfiles = {}
        self.areas = {}
        self.exposures = {}
        self.catnames = {}
        
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
            print(f'SExtractor initialized\nField: {self.field}\nZeropoint: {self.zeropoint}\nConfig File:' +
                  f'{self.config_file}\nCatalog Directory: {self.cat_dir}\nImage Directory: {self.image_dir}' +
                  f'\nImage File Pattern: {self.imfile}')

    def auto_psf(self): #default psf dictionaries
        self.psf_fwhm = {'F090W': 0.06, 'F115W': 0.05, 'F150W': 0.055, 'F182M': 0.06, 'F200W': 0.065, 
                         'F210M': 0.07, 'F277W': 0.1, 'F300M': 0.11, 'F335M': 0.12, 'F356W': 0.15, 
                         'F360M': 0.15, 'F410M': 0.15, 'F444W': 0.15}
        
        self.psf_max = {'F090W': 0.11, 'F115W': 0.11, 'F150W': 0.08, 'F182M': 0.085, 'F200W': 0.09,
                        'F210M': 0.1, 'F277W': 0.13, 'F300M': 0.14, 'F335M': 0.15, 'F356W': 0.18,
                        'F360M': 0.17, 'F410M': 0.17, 'F444W': 0.17}
        
        self.psf_mag = {'F090W': 26, 'F115W': 26, 'F150W': 26, 'F182M': 26, 'F200W': 26,
                        'F210M': 26, 'F277W': 26, 'F300M': 26, 'F335M': 26, 'F356W': 26,
                        'F360M': 26, 'F410M': 26, 'F444W': 26}
    
    def auto_sb(self): #default sensitivity dictionaries
        self.sb_min = {'F090W': 29.6, 'F115W': 29.4, 'F150W': 29.5, 'F182M': 28.1, 'F200W': 29.7,
                       'F210M': 27.8, 'F277W': 30.3, 'F300M': 29, 'F335M': 29, 'F356W': 30.3,
                       'F360M': 28.9, 'F410M': 30, 'F444W': 30}
        
        self.sb_turn = {'F090W': 0.2, 'F115W': 0.18, 'F150W': 0.21, 'F182M': 0.16, 'F200W': 0.19,
                        'F210M': 0.17, 'F277W': 0.21, 'F300M': 0.18, 'F335M': 0.18, 'F356W': 0.22,
                        'F360M': 0.18, 'F410M': 0.2, 'F444W': 0.2}
        
        self.sb_slope = {'F090W': -5, 'F115W': -4.8, 'F150W': -5, 'F182M': -5, 'F200W': -5.1,
                         'F210M': -5, 'F277W': -5, 'F300M': -5, 'F335M': -5, 'F356W': -4.9,
                         'F360M': -5, 'F410M': -4.9, 'F444W': -5}
            
    def set_psf(self, filt, fwhm, end, mag): #minimum and maximum FWHM to classify an object as a star
        self.psf_fwhm[filt] = fwhm
        self.psf_max[filt] = end
        self.psf_mag[filt] = mag
        
    def set_sb(self, filt, minimum, turn, slope): #sets surface brightness / point source sensitivity curve
        self.sb_min[filt] = minimum
        self.sb_turn[filt] = turn
        self.sb_slope[filt] = slope
        
    def get_images(self): #get an invdividual image and weight file
        #imfile should be a filename pattern (e.g., 'jwst_field_*_30mas_sci.fits')
        #the filter should be replaced by * to get all files of that pattern
        #the file pattern should have a 'wht' file for each 'sci' (or 'drz') file
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
        
    def get_pixel(self):
        for i in range(len(self.imfiles)):
            file = list(self.imfiles.values())[i]
            filt = list(self.imfiles.keys())[i]
            hdul = fits.open(file)
            header = hdul[0].header
                        
            if int(filt[1:4]) < 250: #related to central wavelength (SW)
                self.sw_pixel = header['PIXAR_A2']**0.5
                
            else:
                self.lw_pixel = header['PIXAR_A2']**0.5
                
            hdul.close()
            
        print(f'Pixel sizes: SW-{self.sw_pixel} LW-{self.lw_pixel}')
    
    def get_areas(self): #used for number counts    
        for i in range(len(self.imfiles)):
            file = list(self.imfiles.values())[i]
            filt = list(self.imfiles.keys())[i]
            hdul = fits.open(file)
            data = hdul[0].data
            
            idx = np.where(data > 0)
            
            if int(filt[1:4]) < 250:
                self.areas[filt] = len(idx[0]) * self.sw_pixel**2 / 60**4
            else:
                self.areas[filt] = len(idx[0]) * self.lw_pixel**2 / 60**4
            
            hdul.close()
            
        print('Areas calculated.')
        
    def get_exposures(self): #input to source extractor as gain = exposure time
        for i in range(len(self.imfiles)):
            file = list(self.imfiles.values())[i]
            filt = list(self.imfiles.keys())[i]
            hdul = fits.open(file)
            header = hdul[0].header
                        
            if int(filt[1:4]) < 250: #related to central wavelength (SW)
                self.exposures[filt] = header['XPOSURE'] / 8
                
            elif int(filt[1:4]) >= 250: #related to central wavelength (LW)
                self.exposures[filt] = header['XPOSURE'] / 2
                
            else:
                print('header matching failed')
                
            hdul.close()
            
        print('Exposure times found.')
    
    def sextract(self, dual=False, overwrite=True): #run source extractor; dual=False for single image
        for i in range(len(self.imfiles)): #dual=valid filter for dual image (filter to use as detection image)
            catname = os.path.join(self.cat_dir, f'{self.field}_{list(self.imfiles.keys())[i]}_cat.txt')
            filt = list(self.imfiles.keys())[i]
            self.catnames[filt] = catname
            
            if os.path.exists(catname) and overwrite == True:
                os.system(f'rm {catname}')
            
            if not os.path.exists(catname):
                if dual == False:
                    os.system(f'sex {self.imfiles[filt]} -c {self.config_file} -MAG_ZEROPOINT ' +
                              f'{self.zeropoint} -CATALOG_NAME {catname} -GAIN {self.exposures[filt]} ' +
                              f'-WEIGHT_IMAGE {self.whtfiles[filt]} -WEIGHT_TYPE MAP_WEIGHT')
                else:
                    os.system(f'sex {self.imfiles[dual]} {self.imfiles[filt]} -c {self.config_file} -MAG_ZEROPOINT ' +
                              f'{self.zeropoint} -CATALOG_NAME {catname} -GAIN {self.exposures[filt]} ' +
                              f'-WEIGHT_IMAGE {self.whtfiles[filt]} -WEIGHT_TYPE MAP_WEIGHT')
                
        print('Source Extractor done.')
        
    def go(self, dual=False, ow=False): #run code after running self.initialize()
        if self.field == 'N/A':
            print('Cannot go without initializing field (call sextractor.initialize())')
            
        self.get_images()
        self.get_pixel()
        self.get_areas()
        self.get_exposures()
        self.sextract(dual=dual, overwrite=ow)
        
    def mags_fwhm(self, table, filt):
        photom = Table()
        
        idx = np.where((table['MAGERR_AUTO'] < 5) & (table['ALPHA_J2000'] != 0) & (table['DELTA_J2000'] != 0))
        
        photom['mag'] = table['MAG_AUTO'][idx]
        photom['err'] = table['MAGERR_AUTO'][idx]
        
        if int(filt[1:4]) < 250:
            photom['fwhm'] = self.sw_pixel * table['FWHM_IMAGE'][idx]
            
        else:
            photom['fwhm'] = self.lw_pixel * table['FWHM_IMAGE'][idx]

        return photom
    
    def plot_counts(self, table, filt):
        area = self.areas[filt]
        cutoff1 = self.psf_fwhm[filt]
        cutoff2 = self.psf_max[filt]
        magcutoff = self.psf_mag[filt]
        
        photom = self.mags_fwhm(table, filt)
        star = Table()
        gal = Table()
        
        idx = np.where((photom['mag'] < magcutoff) & (photom['fwhm'] > cutoff1) & (photom['fwhm'] < cutoff2))
        star['mag'] = photom['mag'][idx]
        star['fwhm'] = photom['fwhm'][idx]
        
        idx = np.where(((photom['mag'] > magcutoff) & (photom['fwhm'] < cutoff2) & (photom['fwhm'] > cutoff1))
                        | (photom['fwhm'] > cutoff2))
        gal['mag'] = photom['mag'][idx]
        gal['fwhm'] = photom['fwhm'][idx]
            
        bins = np.linspace(0, 40, num=81)
        bin_indices1 = list(np.digitize(star['mag'], bins) / 2 - 0.25)
        bin_indices2 = list(np.digitize(gal['mag'], bins) / 2 - 0.25)

        x = np.linspace(0.25, 39.75, num=80)  
        y_vals = Table()
        y_vals['star'] = [bin_indices1.count(i) / area for i in x]
        y_vals['star_err'] = [math.sqrt(bin_indices1.count(i)) / area for i in x]
        y_vals['gal'] = [bin_indices2.count(i) / area for i in x]
        y_vals['gal_err'] = [math.sqrt(bin_indices2.count(i)) / area for i in x]
    
        with open(os.path.join(self.cat_dir, f'{self.field}_{filt}_gal_cts.txt'), 'w+') as f:
            f.write(f'# Area={area} sq deg \n # Col 1 = bin center in AB mag \n # Col 2 = N = Number of galaxes' +
                    f' \n # Col 3 = N/A = number per sq deg \n # Col 4 = sqrt(N)/A = poisson error in N/A')
            
            for i in range(len(x)):
                f.write(f"{x[i]}    {y_vals['gal'][i]*area}    {y_vals['gal'][i]}    {y_vals['gal_err'][i]} \n")
        
        with open(os.path.join(self.cat_dir, f'{self.field}_{filt}_star_cts.txt'), 'w+') as f:
            f.write(f'# Area={area} sq deg \n # Col 1 = bin center in AB mag \n # Col 2 = N = Number of stars ' +
                    f'\n # Col 3 = N/A = number per sq deg \n # Col 4 = sqrt(N)/A = poisson error in N/A')
            
            for i in range(len(x)):
                f.write(f"{x[i]}    {y_vals['star'][i]*area}    {y_vals['star'][i]}    {y_vals['star_err'][i]} \n")

        return x, y_vals, star, gal, photom
    
    def plot_cats(self, table, filt):
        cutoff1 = self.psf_fwhm[filt]
        cutoff2 = self.psf_max[filt]
        magcutoff = self.psf_mag[filt]
        SB_min = self.sb_min[filt]
        SB_turn = self.sb_turn[filt]
        SB_slope = self.sb_slope[filt]
        
        x, y_vals, star, gal, photom = self.plot_counts(table, filt)
    
        x_line = np.linspace(SB_turn, 10, num=100)
        y_line = [SB_min + SB_slope * (math.log10(i) - math.log10(SB_turn)) for i in x_line]
    
        fig, ax = plt.subplots(1,3)    
        fig.set_size_inches(18,5)
        
        ml1 = MultipleLocator(0.05)
        ml2 = MultipleLocator(0.5)
        ml3 = MultipleLocator(1)
        ml4 = MultipleLocator(0.5)
    
        ax[0].scatter(photom['mag'], photom['err'], color ='black', s=0.1, zorder=3)
        ax[0].plot([0, 40], [0, 0], lw=1, color='black', ls='--', zorder=8)
        ax[0].plot([0, 40], [0.2, 0.2], lw=1, color='black', ls='--', zorder=8)
        ax[0].set_title(f'{self.field} / {filt}', fontsize=20)
        ax[0].set_xlabel('AB [Mag]', fontsize=15)
        ax[0].set_ylabel('magerr', fontsize=15)
        ax[0].set_ylim(-0.1, 1)
        ax[0].set_xlim(15, 32)
        ax[0].yaxis.set_minor_locator(ml1)
        ax[0].xaxis.set_minor_locator(ml2)
        ax[0].text(16, 0.22, 'SN~5', c='black', fontsize=12)
    
        ax[1].scatter(star['fwhm'], star['mag'], color = 'red', s=0.1, zorder=4)
        ax[1].scatter(gal['fwhm'], gal['mag'], color = 'blue', s=0.1, zorder=4)
        ax[1].plot([cutoff1, cutoff1], [0, 40], color='black', lw=1)
        ax[1].plot([cutoff2, cutoff2], [0, magcutoff], color='black', lw=1)
        ax[1].plot([cutoff1, cutoff2], [magcutoff, magcutoff], color='black', lw=1)
        ax[1].plot([cutoff1, SB_turn], [SB_min, SB_min], color='green', ls='--', lw=2, zorder=5)
        ax[1].plot(x_line, y_line, color='green', ls='--', lw=2, zorder=5)
        ax[1].set_title(f'{self.field} / {filt}', fontsize=20)
        ax[1].set_xlabel('FWHM [arcsec]', fontsize=15)
        ax[1].set_ylabel('AB [Mag]', fontsize=15)
        ax[1].set_ylim(32, 15)
        ax[1].set_xlim(0.04, 5)
        ax[1].set_xscale('log')
        ax[1].yaxis.set_minor_locator(ml3)
        ax[1].text(1.9, 16, "Stars", c='red', fontsize=12)
        ax[1].text(1.9, 16.7, 'Galaxies', c='blue', fontsize=12)

        ax[2].errorbar(x, y_vals['star'], y_vals['star_err'], color = 'red', fmt='o', capsize=4, ms=5, zorder=3)
        ax[2].errorbar(x, y_vals['gal'], y_vals['gal_err'], color = 'blue', fmt='o', capsize=4, ms=5, zorder=3)
        ax[2].set_title(self.field + ' / ' + filt, fontsize=20)
        ax[2].set_xlabel('AB [Mag]', fontsize=15)
        ax[2].set_ylabel('number/sq deg.', fontsize=15)
        ax[2].set_ylim(10**2, 2*10**6)
        ax[2].set_xlim(15, 32)
        ax[2].set_yscale('log')
        ax[2].xaxis.set_minor_locator(ml4)
        ax[2].text(15.5, 1.1*10**6, "Stars", c='red', fontsize=12)
        ax[2].text(15.5, 7*10**5, 'Galaxies', c='blue', fontsize=12)
    
        for axis in ax:
            axis.tick_params(axis='both', which='major', length=6, width=2, direction='in', labelsize='large', 
                             bottom=True, top=True, left=True, right=True)
            axis.tick_params(axis='both', which='minor', length=4, width=1, direction='in', labelsize='large', 
                             bottom=True, top=True, left=True, right=True)

        fig.savefig(os.path.join(self.cat_dir, f'{self.field}_{filt}_plot.png'))

    def run_plots(self):
        for filt in self.imfiles.keys():
            t = Table.read(self.catnames[filt], format='ascii') 
            
            self.plot_cats(t, filt)

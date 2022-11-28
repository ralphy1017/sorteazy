import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.table import Table
from glob import glob

class Eazy:
    
    ### INITIALIZE

    def __init__(self, eazypath, params={}):
        self.eazypath = eazypath          ### Must be in the inputs directory of eazy-photoz
        self.params = params              ### ex: /Users/rafaelortiz/jwst/eazy-photoz/inputs
    
    ### Method makeparam
    ### uses params and writes them to EAZY parameter file

    def makeparam(self):
        with open(f'{self.eazypath}/zphot.param', 'w') as f:
            for param, value in self.params.items():
                f.write(f'{param} {value}\n')
            f.close()

        print(f'Parameters written to zphot.param at {self.eazypath}/zphot.param')
        
    ### Method run
    ### run eazy
    ### can do it with or without a logfile

    def run(self, logfile=True):
        self.logfile = logfile
        os.chdir(self.eazypath)
        if self.logfile == False:
            os.system('../src/eazy')
        else:
            os.system('../src/eazy > logfile')
        
    ### Method makeplots
    ### main plotting function

    def makePlots(self, output_dir, imrgb, zmin=7, zmax=100):
        self.output_dir = output_dir
        namestyle = 'photz_' # for EAZY output files

        # retreive files and get id numbers
        files = glob(os.path.join(output_dir, namestyle + "*.obs_sed"))
        indices = [i.split('_')[-2].split('.')[0] for i in files]
        
        # get main catalog file as astropy table, filter high-z objects
        z_cat = os.path.join(output_dir, 'photz.zout')
        t = Table.read(z_cat, format='ascii')
        high_z = t[np.where((t['z_p'] > zmin) & (t['z_p'] < zmax))]['id']

        # get SExtractor matched catalog file for object positions
        catfile = self.params['CATALOG_FILE']
        t = Table.read(catfile, format='ascii')
        
        # loop through high-z objects
        for idx in high_z:
            # three EAZY files per object
            obs_sed = os.path.join(output_dir, namestyle+str(idx)+".obs_sed")
            pz = os.path.join(output_dir, namestyle+str(idx)+".pz")
            temp_sed = os.path.join(output_dir, namestyle+str(idx)+".temp_sed")

            # open template sed file
            temp_hl = open(temp_sed, 'r')
            lines = temp_hl.readlines()
            z_a = lines[1].split()[1] # redshift
            z_p = lines[2].split()[1] # redshift with prior

            # sometimes there is a space after z= based on how large the redshift is
            if z_a == 'z=':
                z_a += lines[1].split()[2]
            if z_p == 'z_prior=':
                z_p += lines[2].split()[2]

            temp_hl.close()

            # make seds and P(z) into tables
            obs_sed = Table.read(obs_sed, format='ascii')
            pz = Table.read(pz, format='ascii')
            temp_sed = Table.read(temp_sed, format='ascii')

            # get main parameters from the tables
            lamb_sed = temp_sed['lambda']
            flux_sed = temp_sed['tempflux']
            lamb_img = obs_sed['lambda']
            flux_img = obs_sed['flux_cat']
            ferr_img = obs_sed['err_cat']
            z = pz['z']
            p_z = pz['pz']

            ### filter response widths
            filter_width = []
            for i in lamb_img:
                if i - 4.4e+04 > 0: # F444W
                    filter_width.append(0.553*10**4)
                elif i - 4e+04 > 0: # F410M
                    filter_width.append(0.219*10**4)
                elif i - 3.5e+04 > 0: # F356W
                    filter_width.append(0.42*10**4)
                elif i - 2.7e+04 > 0: # F277W
                    filter_width.append(0.356*10**4)
                elif i - 1.9e+04 > 0: # F200W
                    filter_width.append(0.236*10**4)
                elif i - 1.5e+04 > 0: # F150W
                    filter_width.append(0.169*10**4)
                elif i - 1.15e+04 > 0: # F115W
                    filter_width.append(0.135*10**4)
                elif i - 1.05e+04 > 0: # F090W
                    filter_width.append(0.105*10**4)

            ### Color image
            row = t[np.where(t['id'] == idx)] # sort from table of objects
            xpos = 4465 - int(row['Y']) #rgb image coords are weird
            ypos = int(row['X'])
            rgb = imrgb[(xpos-20):(xpos+20), (ypos-20):(ypos+20)] #make 40x40 pixel cutout
            
            ### Plot
            fig, ax = plt.subplots(1,3)    
            fig.set_size_inches(15,5)
            # minor tick locators
            ml1 = MultipleLocator(0.2*10**4)
            ml2 = MultipleLocator(0.5)

            # first plot (SED)
            # template SED
            ax[0].plot(lamb_sed, flux_sed, lw=1, color='black', ls='-', zorder=8)
            # observed SED (each filter)
            ax[0].errorbar(x=lamb_img, y=flux_img, xerr=filter_width, yerr=ferr_img, color = 'red', fmt='o', capsize=4, ms=5, zorder=3)
            fig.suptitle(f'Galaxy id {idx} EAZY Results | {z_a} {z_p} ', fontsize=20)
            ax[0].set_xlabel('Wavelength [$\AA$]', fontsize=15)
            ax[0].set_ylabel(r'$F_{\nu}$ [nJy]', fontsize=15)
            ax[0].set_xlim(0.5*10**4, 5*10**4)
            ax[0].set_ylim(min(flux_img)*0.2, max(flux_img)*1.5)
            ax[0].xaxis.set_minor_locator(ml1)
            ax[0].set_yscale('log')
            ax[0].set_title('SED Fit', fontsize=15)

            # second plot (probability)
            ax[1].plot(z, p_z, lw=1, color='blue', ls='-')
            ax[1].set_xlabel('z', fontsize=15)
            ax[1].set_ylabel('P',fontsize=15)
            ax[1].xaxis.set_minor_locator(ml2)
            ax[1].set_title('P(z)', fontsize=15)

            # third plot (color image)
            # just in case there are errors
            try:
                ax[2].imshow(rgb)
            except ValueError:
                ax[2].set_facecolor('black') # just set it black if it fails
            
            # tick parameters
            for axis in ax:
                axis.tick_params(axis='both', which='major', length=6, width=2, direction='in', labelsize='large', 
                                bottom=True, top=True, left=True, right=True)
                axis.tick_params(axis='both', which='minor', length=4, width=1, direction='in', labelsize='large', 
                                bottom=True, top=True, left=True, right=True)
            
            # save
            fig.savefig(os.path.join(output_dir, f'{idx}_EAZY_plot.png'))
 
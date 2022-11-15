import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.table import Table
import math
from glob import glob

class Eazy:
    
    def __init__(self,eazypath, params={}):
        self.eazypath = eazypath          ### Must be in the inputs directory of eazy-photoz
        self.params = params              ### ex: /Users/rafaelortiz/jwst/eazy-photoz/inputs
    
    def makeparam(self):
        with open(f'{self.eazypath}/zphot.param', 'w') as f:
            for param, value in self.params.items():
                f.write(f'{param} {value}\n')
            f.close()
        
    def run(self):
        os.chdir(self.eazypath)
        os.system('../src/eazy > logfile')
        
    def makePlots(self, output_dir, namestyle):
        self.output_dir = output_dir
        self.namestyle = namestyle

        files = glob(os.path.join(output_dir, namestyle + "*.obs_sed"))
        indices = [i.split('_')[-2].split('.')[0] for i in files]
        idx = indices[39]
        
        z_cat = os.path.join(output_dir, 'photz.zout')
        t = Table.read(z_cat, format='ascii')
        high_z = t[np.where(t['z_p'] > 5)]['id']
        odds = {}
        for i in range(len(t['id'])):
            odds[i+1] = t[i]['odds']
        
        for idx in high_z:
            obs_sed = os.path.join(output_dir, namestyle+str(idx)+".obs_sed")
            pz = os.path.join(output_dir, namestyle+str(idx)+".pz")
            temp_sed = os.path.join(output_dir, namestyle+str(idx)+".temp_sed")
            temp_hl = open(temp_sed, 'r')
            lines = temp_hl.readlines()
            z_a = lines[1].split()[1]
            z_p = lines[2].split()[1]
            if z_a == 'z=':
                z_a += lines[1].split()[2]
            if z_p == 'z_prior=':
                z_p += lines[2].split()[2]
            temp_hl.close()
            obs_sed = Table.read(obs_sed, format='ascii')
            pz = Table.read(pz, format='ascii')
            temp_sed = Table.read(temp_sed, format='ascii')

            lamb_sed = temp_sed['lambda']
            flux_sed = temp_sed['tempflux']
            lamb_img = obs_sed['lambda']
            flux_img = obs_sed['flux_cat']
            ferr_img = obs_sed['err_cat']
            z = pz['z']
            chi2 = pz['pz']
            odd = odds[int(i)]
            
            fig, ax = plt.subplots(1,2)    
            fig.set_size_inches(12,6)
            ml1 = MultipleLocator(0.2*10**4)
            ml3 = MultipleLocator(0.5)
            
            ax[0].plot(lamb_sed, flux_sed, lw=1, color='black', ls='-', zorder=8)
            ax[0].errorbar(x=lamb_img, y=flux_img, yerr=ferr_img, color = 'red', fmt='o', capsize=4, ms=5, zorder=3)
            fig.suptitle(f'Star id {idx} EAZY Results | {z_a} {z_p} | Q$_z$ = {odd}', fontsize=20)
            ax[0].set_xlabel('Wavelength [Angstroms]', fontsize=15)
            ax[0].set_ylabel('Flux [F$_{\lambda}$]', fontsize=15)
            ax[0].set_xlim(0.5*10**4, 4.7*10**4)
            ax[0].set_ylim(min(flux_img)*0.4 if min(flux_img) > 5 else 0, max(flux_img)*1.5)
            ax[0].xaxis.set_minor_locator(ml1)
            ax[0].set_title('SED Fit', fontsize=15)
            ax[1].plot(z, chi2, lw=1, color='blue', ls='-')
            ax[1].set_xlabel('z', fontsize=15)
            ax[1].set_ylabel('P',fontsize=15)
            ax[1].xaxis.set_minor_locator(ml3)
            ax[1].set_title('P(z)', fontsize=15)
            
            for axis in ax:
                axis.tick_params(axis='both', which='major', length=6, width=2, direction='in', labelsize='large', 
                                bottom=True, top=True, left=True, right=True)
                axis.tick_params(axis='both', which='minor', length=4, width=1, direction='in', labelsize='large', 
                                bottom=True, top=True, left=True, right=True)
            
            fig.savefig(os.path.join(output_dir, f'{idx}_EAZY_plot.png'))
 


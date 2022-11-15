from sextractor import Sextractor
import os
from glob import glob

elgordo = Sextractor("elgordo",28.0865,'/Users/rafaelortiz/sextractor/config.sex', 
                        '/Users/rafaelortiz/Work/jwst/pearls/el-gordo/30mas/catalogs', 
                        '/Users/rafaelortiz/Work/jwst/pearls/el-gordo/30mas', 
                        'mosaic_elgordo_total_nircam_*_30mas_20221014_drz.fits',
                        dual=False, overwrite=True)

elgordo.run_plots('/Users/rafaelortiz/Work/jwst/pearls/el-gordo/30mas/catalogs')
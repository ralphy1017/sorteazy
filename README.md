# sorteazy
Perform photometry on JWST images using Source Extractor and EAZY.

## Important Classes and Methods

### prepimg

Convert default MAST i2d.fits files into pixel-aligned science and weight extensions.

__init__() - Splits and reprojects i2d.fits files and prints the field, image directory, and image naming scheme, which can be fed to the sextractor class.

           Arguments: 
           image_dir - directory in which the i2d.fits files are stored
           file - pattern of filenames of the images (replace the filter with *)        
           field - field name to put in image and catalog files

           Returns: Nothing, but prints field, image_dir, new file pattern, and zeropoints


bkgsub() - Can be called after go() to subtract the background from the images.

           Arguments: 
           size - size of tiles to compute and subtract background from (default 100 pixels)

           Returns: Nothing, but prints field, image_dir, and new file pattern


### sextractor

Run source extractor on science and weight images created by prepimg, and produce plots.

__init__() - gets necessary info to run source extractor

           Arguments: 
           field - field name
           image_dir - directory that contains images
           imfile - image file naming scheme
           cat_dir - directory to place catalogs
           zp - array of zeropoints in the form (zp_sw, zp_lw) (default 28.0865 for both - true for 30mas pixels)

           Returns: Nothing

set_psf() - sets bounds of fwhm and magnitude to do star-galaxy separation

           Arguments: 
           filt - filter that it applies to
           fwhm - minimum fwhm that an object (star or galaxy) can have to be counted in the plots
           end - maximum fwhm that a star can have
           mag - maximum magnitude that a star can have

           Returns: Nothing

set_sb() - sets shape of surface brightness and point source sensitivity limit

           Arguments: 
           filt - filter that it applies to
           sb_min - point source sensitivity limit (magnitude)
           sb_turn - fwhm where the point source and surface brightness sensitivity limits meet
           sb_slope - slope of the surface brightness sensitivity

           Returns: Nothing

sextract() - gets images, gets pixel sizes, gets areas, gets exposure times, and runs source extractor

           Arguments: 
           config_file - absolute path to SExtractor configuration file
           dual - dual image mode option (default False=don't run dual image, else specify filter to use as detection image)
           overwrite - whether or not to overwrite the files (default False)

           Returns: Nothing

run_plots() - creates photometry plots using previously created SExtractor catalogs

           Arguments: None
           Returns: Nothing
           

# sorteazy
Perform photometry on JWST images using Source Extractor and EAZY.

## Important Classes and Methods

### PrepImg.py

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


### SExtractor.py

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
           

### MakeCats.py

Makes a custom catalog in AB Mag using each catalog generated from SExtractor so that it may be passed to eazy to find the Photometric Redshift of each object in the catalog

__init__() - Receives necessary info to initialize a field with its catalogs
            
            Arguments:
            field - the field name
            cat_dir - catalog directory where SExtractor catalogs exist on your machine
            cat_file - catalog file namestyle so that it can group catalogs in each filter
            pixsize - pixel size found in fits file data

	    Returns: Nothing

match() - Filters nondetections and aligns each catalog using SkyCoord and then writes to a single catalog for eazy

            Arguments:
            det_filt - string of the name of the filter to use as the matching filter

	    Returns: Nothing


### eazyMethods.py

Runs eazy wherever it exists on your machine and then generates SED vs Chi^2 Plots with eazy photometric Output

__init__() - Receives location of eazy and the custom parameters 

            Arguments:
            eazypath - path to eazy-photoz on your machine to navigate files properly
            params - dictionary of parameters with correct parameter formatting to then write to zphot.param so that it  
                     overides zphot.param.default

	    Returns: Nothing

makeparam() - makes zphot.param file 

	    Arguments: None
	    Returns: Nothing

run() - runs eazy on your machine, set logfile=False to show output. Otherwise, output is written to a log file 
        ex: ../src/eazy > logifle

	    Arguments: 
	    logfile - logfile=False to show output, logfile=True to save output to a log file (default True)

	    Returns: Nothing (but runs EAZY so saves all of the resulting files in OUTPUT directory)

makePlots() - makes SED template plots, Chi^2 vs Z plots from eazy output, and color stamps of individual high-z galaxies. 

            Arguments:
            output_dir - path of eazy OUTPUT directory
	    imrgb - array of RGB values for entire field
	    zmin - minimum redshift value to make plots for (default 7)

	    Returns: eazy plots + color stamps of high-z galaxies


### Trilogy.py

Creates color images from science images so that it can be fed to eazyMethods.py for plotting individual galaxies.

__init__() - receives location of image files, as well as the field name for saving images.

	    Arguments:
	    image_dir - directory that contains science images
	    imfiles - pattern of filenames within image_dir
	    field - name of field

	    Returns: Nothing, but prints filters and images corresponding to each filter

bkgsub() - subtracts background from science images

	    Arguments: None
	    Returns: Nothing

auto_colors() - assigns colors to each filter from purple (short wavelength) to red (long wavelength)

	    Arguments: None
	    Returns: prints contribution of each filter to R, G, B (in that order)

make_stamp() - creates a color image stamp to assign levels to the various parameters of Trilogy.

	    Arguments:
	    sample_size - side width of square image (default 1000)
	    dx - shift from the center of the image in pixels (right = positive) (default 0)
	    dy - shift from the center of the image in pixels (up = positive) (default 0)
	    noiselum - saturation value of pixels at noisesig sigma above the mean surface brightness value (default 0.12) (spans from 0-1)
	    satpercent - percent of pixels to be saturated at 1.0 (default 0.001)
 	    noisesig - number of sigma above the mean for pixels to have a value of noiselum (default 1)
	    noisefloorsig - number of sigma below the mean for pixels to be black (0.0) (default 2)
	    
	    Returns: color image stamp as a matplotlib plot

make_RGB() - makes full field RGB image from scaling found with make_stamp()

	    Arguments: None
	    Returns: Nothing

enhance() - enhances full field RGB image and saves it.

	    Arguments:
	    color - color enhancement multiplier (1.0 = no change) (default 1.3)
	    contrast - contrast enhancement multiplier (default 1.4)
	    brightness - brightness enhancement multiplier (default 1.1)
	    sharpness - sharpness enhancement multiplier (default 1.9)

	    Returns: array of RGB values, and saves enhanced image in image_dir
		    

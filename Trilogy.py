# basic imports
import numpy as np
import os
from glob import glob
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
# to suppress warnings
import warnings

# astropy and photutils and image things
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import photutils  # for background subtraction
from scipy.optimize import golden
from PIL import Image, ImageEnhance

class Trilogy:

    ### INITIALIZE

    def __init__(self, image_dir, imfiles, field):
        # ignore warnings
        warnings.simplefilter('ignore')

        # make empty data structures
        self.images = {}
        self.filters = []
        self.image_data = {}
        self.bkgsub_data = {}

        # get image dir and image file pattern
        self.image_dir = image_dir
        self.imfiles = imfiles
        self.field = field

        # retreive images
        self.imfiles = glob(os.path.join(image_dir, imfiles))
        self.imfiles.sort()

        # get filters and make dictionary of image files
        for image in self.imfiles:
            imfile = image.lower().split('/')[-1]
            
            for i in range(len(imfile)):
                if imfile[i] == 'f' and (imfile[i+4] == 'w' or imfile[i+4] == 'm') and imfile[i+1] in '0 1 2 3 4 5 6 7 8 9'.split():
                    filt = imfile[i:i+5]
                    break

            self.filters.append(filt.upper())
            self.images[filt.upper()] = image

        # get image data and put in dictionary
        for filt in self.filters:
            self.image_data[filt] = fits.open(self.images[filt])['sci'].data

        # print
        for k, v in zip(self.filters, self.images.values()):
            print(k, v)

    #########################################################################################################
    # START OF TRILOGY METHODS
    #########################################################################################################
    # the only change from Dan Coe's code is putting 'self' before called methods

    ### TRILOGY METHOD 1 (from Dan Coe)

    def da(self, k):
        a1 = k * (x1 - x0) + 1
        a2 = k * (x2 - x0) + 1
        a1n = a1**n
        a1n = np.abs(a1n)  # Don't want the solutions where a1 & a2 are both negative!
        da1 = a1n - a2
        k = np.abs(k)
        if k == 0:
            return self.da(1e-10)
        else:
            da1 = da1 / k  # To avoid solution k = 0!
        return abs(da1)

    ### TRILOGY METHOD 2 (from Dan Coe)

    def imscale2(self, data, levels, y1):
        # x0, x1, x2  YIELD  0, y1, 1,  RESPECTIVELY
        # y1 = noiselum
        global n, x0, x1, x2  # So that golden can use them
        x0, x1, x2 = levels  
        if y1 == 0.5:
            k = (x2 - 2 * x1 + x0) / float(x1 - x0) ** 2
        else:
            n = 1 / y1
            k = np.abs(golden(self.da))
        r1 = np.log10( k * (x2 - x0) + 1)
        v = np.ravel(data)
        v = self.clip2(v, 0, None)
        d = k * (v - x0) + 1
        d = self.clip2(d, 1e-30, None)
        z = np.log10(d) / r1
        z = np.clip(z, 0, 1)
        z.shape = data.shape
        z = z * 255
        z = z.astype(np.uint8)
        return z

    ### TRILOGY METHOD 3 (from Dan Coe)

    def clip2(self, m, m_min=None, m_max=None):
        # nanmin and nanmax important to ignore nan values
        # otherwise you'll get all 0's
        if m_min == None:
            m_min = np.nanmin(m)
        if m_max == None:
            m_max = np.nanmax(m)
        return np.clip(m, m_min, m_max)

    ### TRILOGY METHOD 4 (from Dan Coe)

    def set_levels(self, data, pp, stripneg=False, sortedalready=False):
        if sortedalready:
            vs = data
        else:
            print('sorting...')
            vs = np.sort(data.flat)
        if stripneg:  # Get rid of negative values altogether!
            # This is the way I was doing it for a while
            # Now that I'm not, resulting images should change (get lighter)
            i = np.searchsorted(vs, 0)
            vs = vs[i+1:]
        else:  # Clip negative values to zero
            vs = self.clip2(vs, 0, None)
        ii = np.array(pp) * len(vs)
        ii = ii.astype(int)
        ii = np.clip(ii, 0, len(vs)-1)
        levels = vs.take(ii)
        #print ii, levels, vs, sort(vs)
        return levels

    ### TRILOGY METHOD 5 (from Dan Coe)

    def determine_scaling(self, data, unsatpercent, noisesig=1, correctbias=True, noisefloorsig=2):
        """Determines data values (x0,x1,x2) which will be scaled to (0,noiselum,1)"""
        # Robust mean & standard deviation
        datasorted = data + 0
        datasorted[np.isnan(datasorted)]=0  # set all nan values to zero
        datasorted = np.sort(datasorted.flat)
        if datasorted[0] == datasorted[-1]:  # data is all one value
            levels = 0, 1, 100  # whatever
        else:
            data_mean, data_median, data_stddev = sigma_clipped_stats(datasorted)
            m = data_mean
            r = data_stddev
            #print('%g +/- %g' % (m, r))

            if correctbias:
                x0 = m - noisefloorsig * r
            else:
                x0 = 0
            x1 = m + noisesig * r
            x2 = self.set_levels(datasorted, np.array([unsatpercent]), sortedalready=True)[0]
            levels = x0, x1, x2
        return levels

    ### TRILOGY METHOD 6 (from Dan Coe)

    def stamp_extent(self, data, sample_size=1000, dx=0, dy=0):
        data_shape = data.shape
        if len(data_shape) == 2:
            ny, nx = data.shape
        else:
            ny, nx, three = data.shape    
        #if nx == 16000:
        #    yc = 7482
        #    xc = 7634
        #elif nx == 6400:
        #    yc = 2755
        #    xc = 3417
        yc = int(ny / 2)
        xc = int(nx / 2)
        #print(xc, yc)
        #print(xc+dx, yc+dy)
        
        ylo = yc - sample_size / 2 + dy
        yhi = yc + sample_size / 2 + dy

        xlo = xc - sample_size / 2 + dx
        xhi = xc + sample_size / 2 + dx
        #print(xlo, xhi, ylo, yhi)
        
        ylo = int(np.clip(ylo, 0, ny))
        yhi = int(np.clip(yhi, 0, ny))
        xlo = int(np.clip(xlo, 0, nx))
        xhi = int(np.clip(xhi, 0, nx))
        #print(xlo, xhi, ylo, yhi)
        return xlo, xhi, ylo, yhi

    ### TRILOGY METHOD 7 (from Dan Coe)

    def image_stamps(self, data, sample_size=1000, dx=0, dy=0):
        xlo, xhi, ylo, yhi = self.stamp_extent(data, sample_size, dx, dy)
        stamps = data[ylo:yhi,xlo:xhi]
        return stamps

    #########################################################################################################
    # END OF TRILOGY METHODS
    #########################################################################################################

    ### Method bkgsub
    ### subtracts background from images

    def bkgsub(self, size=100):
        # empty dictionary of background data (arrays)
        background_maps = {}

        # find background maps
        for filt in tqdm(self.filters, desc='Finding background maps'):
            background_maps[filt] = photutils.Background2D(self.image_data[filt], size, filter_size=5) #size - size of bkg boxes

        # subtract background
        for filt in tqdm(self.filters, desc='Subtracting backgrounds'):
            self.bkgsub_data[filt] = self.image_data[filt] - background_maps[filt].background

    ### Method auto_colors
    ### prepares color contribution of each filter to the RGB image

    def auto_colors(self):
        # set colors for many filters
        # it does require boosting saturation after the fact

        cmap = matplotlib.cm.get_cmap("rainbow")

        self.filter_colors = {}
        
        # get filter color contributions
        for i, filt in enumerate(self.filters):
            x_min = 0.0  # bluest filter will be purple
            x = i / (len(self.filters) - 1) * (1 - x_min) + x_min
            r_lum, g_lum, b_lum, alpha = cmap(x)
            rgb_lum = np.array([r_lum, g_lum, b_lum])
            self.filter_colors[filt] = rgb_lum
            print(filt, ' %4.2f' % r_lum, ' %4.2f' % g_lum, ' %4.2f' % b_lum)
            
        ### filter color sum

        # make filter colors into numpy arrays
        for filt in self.filters:
            self.filter_colors[filt] = np.array(self.filter_colors[filt])
            
        # get total contribution to each color (R, G, B)    
        self.rgb_lum_sum = np.zeros(3)
        for i, filt in enumerate(self.filters):
            self.rgb_lum_sum += np.array(self.filter_colors[filt])

    ### Method make_stamp
    ### makes RGB stamp image (and saves how the colors are made)

    def make_stamp(self, sample_size=1000, dx=0, dy=0, noiselum=0.12, satpercent=0.001, noisesig=1, noisefloorsig=2):
        # parameters
        self.sample_size = sample_size
        self.dx = dx
        self.dy = dy
        self.noiselum = noiselum
        self.satpercent = satpercent
        self.unsatpercent = 1 - 0.01 * satpercent
        self.noisesig = noisesig
        self.correctbias = True
        self.noisefloorsig = noisefloorsig

        image_stamp = self.image_stamps

        # make stamp image
        scaled_images = {}
        self.levels_all = {}
        for filt in tqdm(self.filters, desc='Scaling colors'):
            data = self.bkgsub_data[filt]

            my_stamp_extent = self.stamp_extent(data, self.sample_size, self.dx, self.dy)
            stamp = image_stamp(data, self.sample_size, self.dx, self.dy)
            levels = self.determine_scaling(stamp.ravel(), self.unsatpercent, noisesig, self.correctbias, self.noisefloorsig)
            scaled = self.imscale2(stamp, levels, noiselum)
            self.levels_all[filt] = levels
            scaled_images[filt] = scaled
    
        rgb_total = 0
        for filt in tqdm(self.filters, desc='Combining colors'):
            rgb = r, g, b = self.filter_colors[filt][:, np.newaxis, np.newaxis] * scaled_images[filt]
            imrgb = np.array([r, g, b]).transpose((1,2,0)).astype(np.uint8)    
            rgb_total = rgb_total + rgb
        r, g, b = rgb_average = rgb_total / self.rgb_lum_sum[:, np.newaxis, np.newaxis]

        imrgb = np.array([r, g, b]).transpose((1,2,0)).astype(np.uint8)
        fig, ax = plt.subplots(1, 1, figsize=(9.5, 6))
        plt.imshow(imrgb, origin='lower', extent=my_stamp_extent) # (xlo,xhi,ylo,yhi))
    
    ### Method make_RGB
    ### makes RGB image from stamp image scaling

    def make_RGB(self):
        #make image for entire field
        scaled_images = {}
        for filt in tqdm(self.filters, desc='Scaling colors'):
            data = self.bkgsub_data[filt]
            levels = self.levels_all[filt]
            scaled = self.imscale2(data, levels, self.noiselum)
            scaled_images[filt] = scaled
            
        rgb_total = 0
        for filt in tqdm(self.filters, desc='Combining colors'):
            rgb = r, g, b = self.filter_colors[filt][:, np.newaxis, np.newaxis] * scaled_images[filt]
            rgb_total = rgb_total + rgb
            
        r, g, b = rgb_average = rgb_total / self.rgb_lum_sum[:, np.newaxis, np.newaxis]

        self.imrgb = np.array([r, g, b]).transpose((1,2,0)).astype(np.uint8)

        print('Done making full RGB image.')

    ### Method enhance
    ### increases contrast and saturation of image

    def enhance(self, color=1.3, contrast=1.4, brightness=1.1, sharpness=1.9):
        # use PIL
        im = Image.fromarray(self.imrgb, 'RGB')

        # flip image
        im = im.transpose(Image.FLIP_TOP_BOTTOM)

        # enhance
        im = ImageEnhance.Color(im).enhance(color) #1.3 normally
        im = ImageEnhance.Contrast(im).enhance(contrast) #1.4 normally
        im = ImageEnhance.Brightness(im).enhance(brightness) #1.1 normally
        im = ImageEnhance.Sharpness(im).enhance(sharpness) #1.9 normally
        
        # save image
        im.save(os.path.join(self.image_dir, f'{self.field}_RGB_sat.png'))

        self.imrgb_enhanced = np.array(im)
        return self.imrgb_enhanced

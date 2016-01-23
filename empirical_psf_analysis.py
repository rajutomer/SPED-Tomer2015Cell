# empirical_PSF_analysis.py

# This script loads as input eight tif stacks, two for each objective (normal
# and extended).  These tif stacks contain the average bead volume.
# The script outputs eight associated CSV files containing statistics of the
# bead image at each position along the optical axis.

# The script was generated by exporting an ipython notebook.  Original notebook
# cell are demarked below by 'In[#]:'

# Co-author: Aaron Andalman (aaron.andalman@gmail.com)
# Co-author: Vanessa Burns (vburns@gmail.com)

# coding: utf-8


# Import packages and define helper functions
# In[1]:

get_ipython().magic(u'matplotlib inline')
import os, sys, time, datetime, multiprocessing, h5py
from math import ceil, floor
import math
import numpy as np
import scipy as sp
from scipy.ndimage import gaussian_filter, laplace, center_of_mass, label
import scipy.ndimage.filters as filters
from scipy.optimize import curve_fit
import matplotlib.pyplot as pyp
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import tqdm
import pandas as pd

#from IPython.html.widgets import interact

def compute_halfmax_crossings(sig):
    """
    Compute threshold_crossing, linearly interpolated.

    Note this code assumes there is just one peak in the signal.
    """
    half_max = np.max(sig)/2.0
    fwhm_set = np.where(sig > half_max)

    l_ndx = np.min(fwhm_set) #assumes a clean peak.
    if l_ndx > 0:
        fwhm_left_ndx = l_ndx - 1 + ((half_max - sig[l_ndx-1]) / (float(sig[l_ndx]) - sig[l_ndx-1]))
    else:
        fwhm_left_ndx = 0

    r_ndx = np.max(fwhm_set) #assumes a clean peak.
    if r_ndx < len(sig)-1:
        fwhm_right_ndx = r_ndx + ((half_max - sig[r_ndx]) / (float(sig[r_ndx+1]) - sig[r_ndx]))
    else:
        fwhm_right_ndx = len(sig)-1

    return np.array([fwhm_left_ndx,fwhm_right_ndx])

#from http://stackoverflow.com/questions/21242011/most-efficient-way-to-calculate-radial-profile
def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile

def load_image(filename, dtype = None, normalize = False):
    """
    Load the image supplied by the user using OpenCV, and then
    immediately convert it to a numpy array.  The user may request
    that it be cast into a specific data type, or (in the case of
    floating point data) normalized to the range [0, 1.0].
    """
    if not os.path.exists(filename):
        raise IOError("File \"" + filename + "\" does not exist.")

    filetype = filename.split('.')[-1]
    if filetype.lower() == 'tif':
        from libtiff import TIFF
        tif = TIFF.open(filename, mode='r')

        # in order to allocate the numpy array, we must count the directories:
        # code borrowed from TIFF.iter_images():
        depth = 0
        while True:
            depth += 1
            if tif.LastDirectory():
                break
            tif.ReadDirectory()
        tif.SetDirectory(0)

        # Each tiff directory contains one z slice!
        z_count = 0
        for zslice in tif.iter_images():
            # Handle Endian-ness conversion since pylibtiff doesn't do it automatically for some reason.
            if tif.IsByteSwapped():
                zslice = zslice.byteswap()

            # If this is the first slice, allocate the output volume
            if z_count == 0:
                shape = zslice.shape
                im = np.zeros((shape[0], shape[1], depth), dtype=zslice.dtype)

            # Insert this slice
            im[:,:,z_count] = zslice
            z_count += 1

        # If the tiff image has only one dimension, we squeeze it out of existence here.
        if z_count == 1:
            im = np.squeeze(im)

        del tif # Close the image

    elif filetype.lower() == 'jpg':
        # convert RGB to monochromatic
        try:
            import cv2
        except ImportError:
            import cv as cv2
        im = np.asarray(cv2.imread(filename, cv2.CV_LOAD_IMAGE_GRAYSCALE))

    else:
        try:
            import cv2
        except ImportError:
            import cv as cv2
        try:
            im = np.asarray(cv2.imread(filename, -1))
        except:
            im = np.asarray(cv2.LoadImage(filename, -1))
            im = np.asarray(im.ravel()[0][:]) # hack
            print "You are using an old version of openCV. Loading image using cv.LoadImage."

    if not im.shape:
        raise IOError("An error occurred while reading \"" + filename + "\"")

    # The user may specify that the data be returned as a specific
    # type.  Otherwise, it is returned in whatever format it was
    # stored in on disk.
    if dtype:
        im = im.astype(dtype)

    # Normalize the image if requested to do so.  This is only
    # supported for floating point images.
    if normalize :
        if (im.dtype == np.float32 or im.dtype == np.float64):
            return im / im.max()
        else:
            raise NotImplementedError
    else:
        return im


# Specify the datasets to be analyzed
# In[2]:

path = '.'

#Storing metadata as a list. Order is relied upon below.
meta_data = [{'objective':'O4x','medium':'air', 'voxel_size':[0.365, 0.365, 1], 'psffile':'Air_O4x_AvgOf5imagesLo25_Hi100b_rev.Resampled.tif'},
             {'objective':'O4x','medium':'edof', 'voxel_size':[1.46, 1.46, 1], 'psffile':'PSF_O4x_AvgBeadV2_AvgOf10imagesLo40_Hi100b65k_Cent_114Rev_146m1.tif'},
             {'objective':'N10x','medium':'air', 'voxel_size':[0.365, 0.365, 1], 'psffile':'Air_N10x_AvgOf5imagesLo25_Hi100b_rev.Resampled.tif'},
             {'objective':'N10x','medium':'edof', 'voxel_size':[0.65, 0.65, 1], 'psffile':'PSF_N10x_1umBeadsX8avg_25min100max_AffAln_065xy_1z_rev.tif'},
             {'objective':'O10x','medium':'air', 'voxel_size':[0.365, 0.365, 1], 'psffile':'Air_O10X_AvgOf5imagesLo25_Hi100b_rev.Resampled.tif'},
             {'objective':'O10x','medium':'edof', 'voxel_size':[0.585, 0.585, 1], 'psffile':'PSF_O10x_1umBeadsX5avg_25min100max_Aff2Aln_0585xy1z_rev.tif'},
             {'objective':'O20x','medium':'air', 'voxel_size':[0.365, 0.365, 1], 'psffile':'Air_O20x_AvgOf5imagesLo25_Hi100b_rev.Resampled.tif'},
             {'objective':'O20x','medium':'edof', 'voxel_size':[0.2925, 0.2925, 1], 'psffile':'PSF_O20x_1umBeadsX5avg_25min100max_AffAln_02925xy1z_rev.tif'},
            ]

for md in meta_data:
    if not os.path.isfile(os.path.join(path,md['psffile'])):
        print 'Missing file', os.path.join(path,md['psffile'])

print 'Done'


# Load the tif stacks

# In[3]:

psf_stacks = [] #list of stacks with same order as meta_data list
for i, md in enumerate(meta_data):
    psf_stacks.append(load_image(os.path.join(path,md['psffile'])))
print 'Done'


# Upsample the psf_stacks. This can be slow, so output is saved and can be reloaded below assuming selected datasets are unchanged.
# In[ ]:

zoom = np.array([4.0,4.0,1.0]).astype('float') #zoom for each dimension
order = 3 #cubic

import scipy
for i, md in tqdm.tqdm(enumerate(meta_data)):
    psf_stacks[i] = scipy.ndimage.interpolation.zoom(psf_stacks[i], zoom, order=order)
    meta_data[i]['voxel_size'] = np.array(md['voxel_size'])/zoom
    print 'Upsampled volume shape (x,y,z)',  psf_stacks[i].shape
    print 'Upsampled voxel size',  meta_data[i]['voxel_size']


# In[ ]:

np.savez('upsampled_data.npz', psf_stacks=psf_stacks, meta_data=meta_data, zoom=zoom, order=order)
print 'Saved.'


# In[5]:

d = np.load('upsampled_data.npz')
psf_stacks = d['psf_stacks']
meta_data = d['meta_data']
print 'Loaded.'


# Estimate the bead center in each plane of each PSF.
#
# Output is saved and can be reloaded below assuming selected datasets are unchanged.

# In[ ]:

def find_center(vol, iz, z_window):
    zndx = np.clip(np.arange(iz - (z_window/2),iz + (z_window/2)),0,vol.shape[2]-1)
    img = vol[:,:,zndx].mean(axis=2)
    return np.unravel_index(np.argmax(img),img.shape), img

bead_positions = []
for i, md in tqdm.tqdm(enumerate(meta_data)):
    z_window = 50
    bead_positions.append([])
    for iz in range(psf_stacks[i].shape[2]):
        #the bead position is taken as the point of maximum fluorescence after
        #smoothing the volume along the optical axis using a boxcar filter.
        bp, _ = find_center(psf_stacks[i],iz,z_window=z_window)
        bead_positions[i].append(bp)
    bead_positions[i] = np.array(bead_positions[i])
np.savez(os.path.join(path,'beadpositions.npz'), bead_positions=bead_positions)
print 'Saved.'


# In[7]:

ld = np.load(os.path.join(path,'beadpositions.npz'))
bead_positions = ld['bead_positions']
print 'Loaded bead_positions.'


# Define the range along the optical axis within which the bead is reliably detectable.

# In[8]:

meta_data[0]['valid_zndx'] = np.arange(50,200) #
meta_data[1]['valid_zndx'] = np.arange(50,800) #50,800
meta_data[2]['valid_zndx'] = np.arange(113,200)
meta_data[3]['valid_zndx'] = np.arange(50,535)
meta_data[4]['valid_zndx'] = np.arange(110,170)
meta_data[5]['valid_zndx'] = np.arange(10,483)
meta_data[6]['valid_zndx'] = np.arange(100,170)
meta_data[7]['valid_zndx'] = np.arange(10,409)

clean_range = [[60,150],[40,790],[120,160],[50,500],[120,155],[10,450],[140,165],[25,400]]


# Define function to measure fwhm and various other related statistics.

# In[9]:

def get_raw_fwhm_info(img, center, voxel_size):
    """
    Measure the fwhm of the bead in particular slice (img).

    Note, this code assumes voxels are square in xy-plane.
    """

    win = 1

    #Measure FWHM using orthogonal cross sections.
    cross_sect_x = img[center[0]-win:center[0]+win,:].mean(axis=0)
    halfmax_range_x = compute_halfmax_crossings(cross_sect_x)
    fwhm_cross_sect_x =  np.diff(halfmax_range_x) * voxel_size[0]

    cross_sect_y = img[:,center[1]-win:center[1]+win].mean(axis=1)
    halfmax_range_y = compute_halfmax_crossings(cross_sect_y)
    fwhm_cross_sect_y =  np.diff(halfmax_range_y) * voxel_size[0]

    #Measure FWHM using a crude radial average (not used below).
    radial_avg = radial_profile(img, center=center)
    fwhm_radial = 2* np.diff(compute_halfmax_crossings(radial_avg)) * voxel_size[0]

    #Compute some additional stats on the slice
    img = np.copy(img)
    img[img<0] = 0
    img_peak_energy = img[center[0],center[1]]
    img_total_energy = np.sum(img)
    img_fwhm_energy = np.sum(img[img>(img_peak_energy/2.0)])

    #Compute some additional stats on the cross-section
    cross_total_energy = (np.sum(np.clip(cross_sect_x,0,np.inf)) + np.sum(np.clip(cross_sect_y,0,np.inf)))/2.0
    cross_peak_energy = img[center[0],center[1]]
    def fwhm_energy(signal, h):
        return (np.sum(signal[ceil(h[0]):floor(h[1])]) +
                signal[floor(h[0])] * (ceil(h[0]) - h[0]) +
                signal[ceil(h[0])] * (h[1] - floor(h[1])))
    cross_fwhm_energy = (fwhm_energy(cross_sect_x, halfmax_range_x) +
                         fwhm_energy(cross_sect_y, halfmax_range_y)) / 2.0

    return (fwhm_cross_sect_x, fwhm_cross_sect_y, fwhm_radial,
            cross_sect_x, cross_sect_y, radial_avg,
            img_total_energy, img_peak_energy, img_fwhm_energy,
            cross_total_energy, cross_peak_energy, cross_fwhm_energy)

print 'Done.'


# Estimate the background level of each dataset by examining a region far from the bead.

# In[10]:

raw_fwhms = []
for i,md in tqdm.tqdm(enumerate(meta_data)):
    raw_fwhms.append(pd.DataFrame(columns = ['z_ndx','z_um','max','fwhm_cross_sect','fwhm_radial','central_energy','total_energy','cross_peak_energy'] ))
    for df_ndx,zndx in enumerate(md['valid_zndx']):

        (fwhm_cross_sect_x, fwhm_cross_sect_y, fwhm_radial,
        cross_sect_x, cross_sect_y,radial_avg,total,peak,_,_,cross_peak_energy,_) = get_raw_fwhm_info(psf_stacks[i][:,:,zndx], bead_positions[i][zndx,:], md['voxel_size'])
        fwhm_cross_sect = np.mean([fwhm_cross_sect_x, fwhm_cross_sect_y])
        raw_fwhms[i].loc[df_ndx] = [zndx, zndx*md['voxel_size'][2],radial_avg.max(),fwhm_cross_sect,fwhm_radial,peak,total,cross_peak_energy]

background_level = np.zeros(len(meta_data))
background_std = np.zeros(len(meta_data))
for i, md in enumerate(meta_data):
    s_z = clean_range[i][0]
    e_z = clean_range[i][1]
    f = raw_fwhms[i].query('z_ndx > @s_z and z_ndx < @e_z')
    #min_ndx = np.argmin(f['fwhm_cross_sect'])
    min_ndx = np.argmax(f['cross_peak_energy'])
    zndx = f['z_ndx'][min_ndx]
    img = psf_stacks[i][:,:,zndx]
    sh = img.shape
    if i != 4:
        img = np.copy(img[20:sh[0]/8,sh[1]/2.0 - sh[1]/8.0:sh[1]/2.0 + sh[1]/8.0])
    else:
        #one of the beads is cropped differently, so use a different region to estimate noise.
        img = np.copy(img[bead_positions[i][zndx,0]-sh[0]/4.0:bead_positions[i][zndx,0]+sh[0]/4.0,
                          bead_positions[i][zndx,1]-sh[1]/4.0:bead_positions[i][zndx,1]+sh[1]/4.0])
        sh = img.shape
        img[sh[0]/2.0 - sh[0]/4.0:sh[0]/2.0 + sh[0]/4.0,
            sh[1]/2.0 - sh[1]/4.0:sh[1]/2.0 + sh[1]/4.0] = np.nan
    background_level[i] = np.nanmean(img)
    background_std[i] = np.nanstd(img)
print 'Done.'


# Measure FWHM and other statistics for every slice in every dataset.

# In[11]:

crosses = []
fwhms = []
for i,md in tqdm.tqdm(enumerate(meta_data)):
    fwhms.append(pd.DataFrame(columns = ['z_ndx','z_um','fwhm_cross_sect'] ))
    crosses.append([])
    for df_ndx,zndx in enumerate(md['valid_zndx']):
        #compute fwhm from x cross section of bead
        (fwhm_cross_sect_x, fwhm_cross_sect_y, fwhm_radial,
        cross_sect_x, cross_sect_y,radial_avg,
        img_total,img_peak,img_fwhm,
        cross_total,cross_peak,cross_fwhm) = get_raw_fwhm_info(psf_stacks[i][:,:,zndx] - background_level[i], bead_positions[i][zndx,:], md['voxel_size'])
        fwhm_cross_sect = np.mean([fwhm_cross_sect_x, fwhm_cross_sect_y])
        row = [zndx, zndx*md['voxel_size'][2],fwhm_cross_sect]
        fwhms[i].loc[df_ndx] = [zndx, zndx*md['voxel_size'][2],fwhm_cross_sect]
        crosses[i].append(cross_sect_x)
    crosses[i] = np.array(crosses[i])
print 'Done.'


# Save statistics to csv files for plotting.

# In[13]:

for i, f in enumerate(fwhms):
    f.to_csv('psfdata_%s_%s.csv'%(meta_data[i]['objective'], meta_data[i]['medium']))
print 'Saved.'
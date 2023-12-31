# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:00:12 2020

@author: soren
"""

import numpy as np
import matplotlib.pyplot as plt
import imageio

def readimages(filenameBase):
  ''' Read in images for flat field correction. '''
  # filenameBase: basename of files, to which are added (_dark or _light).tiff
  dark = imageio.imread(filenameBase + '_dark.tiff')
  light = imageio.imread(filenameBase + '_light.tiff')
  return light - dark, dark

def fitprofile(img, h):
  ''' Fit a Gaussian to img, subsampling by h. '''
  s = img.shape
  x = np.arange(0, s[1], h)
  y = np.arange(0, s[0], h)
  data = img[0:s[0]:h, 0:s[1]:h].flatten(order='F')
  data[~np.isfinite(data)] = 0
  ind = (data > 0)
  xx = np.expand_dims(np.repeat(x, len(y)), axis=1)
  yy = np.expand_dims(np.tile(y, len(x)), axis=1)
  xx = xx[ind]
  yy = yy[ind]
  A = np.hstack((xx**2+yy**2, xx, yy, np.ones((len(xx), 1))))
  p = np.linalg.pinv(A) @ np.log(data[ind])
  a = -p[0]
  x0 = p[1]/(2*a)
  y0 = p[2]/(2*a)
  c = np.exp(p[3]+a*(x0**2+y0**2))

  return a, x0, y0, c

def getgaussian(img, r, h, d, display):
  ''' Fit image to a 2D Gaussian of the form exp(-a((x-x0)**2 + (y-y0)**2)) with mean=1
      and return parameters and image of Gaussian. '''
  # Input:
  # img: image to be fit
  # r: radius of pixels to be fit
  # h: spacing of grid of pixels used for fit
  # d: mechanical tolerance in pixels for center of light distribution on sensor
  # display: display figure if True
  # Output:
  # out: image of Gaussian
  # a: Gaussian width parameter
  # x0, y0: Gaussian center position
  # c: Gaussian amplitude

  # approximate center of distribution
  s = img.shape
  y0 = int(s[0]/2)
  x0 = int(s[1]/2)
  mesh = np.mgrid[0:s[0], 0:s[1]]
  img[((mesh[0, :, :]-y0)**2 + (mesh[1, :, :]-x0)**2)>((r+d)**2)] = np.nan
  a, x0, y0, c = fitprofile(img, h)

  # get final Gaussian parameters
  img[((mesh[0, :, :]-y0)**2 + (mesh[1, :, :]-x0)**2)>(r**2)] = np.nan
  a, x0, y0, c = fitprofile(img, h)

  # make Gaussian image with mask
  cal = c*np.exp( -a*( (mesh[0, :, :]-y0)**2 + (mesh[1, :, :]-x0)**2 ))
  cal[((mesh[0, :, :]-y0)**2 + (mesh[1, :, :]-x0)**2)>r**2] = np.nan
  m = np.nanmean(cal)
  cal /= m  # normalize so that mean value = 1

  if display:
    displayfigs(img, m*cal)

  return cal, a, x0, y0, c

def displayfigs(img, cal):
  ''' Display results of calibration. '''
  s = img.shape

  _, ax = plt.subplots(nrows=2, ncols=2)
  ax[0, 0].imshow(img)
  ax[0, 0].set_xticks([])
  ax[0, 0].set_yticks([])
  ax[0, 0].set_title('Original Image')
  ax[0, 1].imshow(cal)
  ax[0, 1].set_xticks([])
  ax[0, 1].set_yticks([])
  ax[0, 1].set_title('Gaussian')
  ax[1, 0].imshow(img-cal)
  ax[1, 0].set_xticks([])
  ax[1, 0].set_yticks([])
  ax[1, 0].set_title('Residual')
  ax[1, 1].imshow(img/cal)
  ax[1, 1].set_xticks([])
  ax[1, 1].set_yticks([])
  ax[1, 1].set_title('Flattened Image')

  _, ax = plt.subplots(nrows=2, ncols=2)
  ax[0, 0].plot(np.arange(0, s[0]), img[:, int(s[1]/2)])
  ax[0, 0].plot(np.arange(0, s[0]), cal[:, int(s[1]/2)])
  ax[0, 0].set_title('X Profile (Original)')
  ax[1, 0].plot(np.arange(0, s[1]), img[int(s[0]/2), :])
  ax[1, 0].plot(np.arange(0, s[1]), cal[int(s[0]/2), :])
  ax[1, 0].set_title('Y Profile (Original)')
  ax[0, 1].plot(np.arange(0, s[0]), np.nanmean(cal)*img[:, int(s[1]/2)]/cal[:, int(s[1]/2)])
  ax[0, 1].set_title('X Profile (Flat)')
  ax[1, 1].plot(np.arange(0, s[1]), np.nanmean(cal)*img[int(s[0]/2), :]/cal[int(s[0]/2), :])
  ax[1, 1].set_title('Y Profile (Flat)')

def hotpixels(dark, nsigmas):
  ''' Generate mask with locations of pixels that have high values at low light levels. '''
  # dark: average of many dark images
  # nsigmas: threshold for defining bad pixels is nsigma standard deviations above mean
  # mask: images of bad pixels

  cut = np.nanmean(dark) + nsigmas*np.nanstd(dark)
  mask = np.ones(dark.shape)
  mask[dark > cut] = np.nan

  return mask

def flattenimg(img, dark, cal, mask):
  ''' Return flattened image with hot pixels removed. '''
  # img: input image
  # dark: dark image
  # cal: image for flat fielding
  # mask: mask of hot pixels

  return (mask/cal)*(img - np.nanmean(mask*dark))

def getstats(img, dark, K):
  ''' Calculate mean, standard deviation, and contrast. '''
  # img: input image (flat w/ dark mean subtracted & hot pixels removed)
  # dark: dark image (hot pixels removed)
  mean = np.nanmean(img)
  var = np.nanvar(img) - np.nanvar(dark) - K*mean
  if var < 0:
    sigma = np.nan #this should not happen
  else:
    sigma = np.sqrt(var)
  contrast = sigma / mean

  return sigma, mean, contrast

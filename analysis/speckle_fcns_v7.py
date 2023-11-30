# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 08:00:15 2020

@author: soren
"""
import json
import csv
import numpy as np
import re
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display
import copy
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
from scipy.ndimage import gaussian_filter1d


def readlogfile(fname, nlocations):

    my_file = open(fname, "r")
    content = my_file.read()
    my_file.close()
    content_list = content.split("\n")

    location = 0
    goodscan = [True]*nlocations
    nattempts = [0]*nlocations

    for i in range(len(content_list)):
        if re.search('Capture at point', content_list[i]):
            #print(content_list[i])
            location +=1
        if re.search('Capture complete', content_list[i]):
            #print(content_list[i])
            nattempts[location-1] +=1
        if re.search('WARNING: Location', content_list[i]):
            print(content_list[i])
            goodscan[location-1] = False

    return goodscan, nattempts



def read_bloodflow_protocsv2(fname):

    # initialize variables
    with open(fname, mode='r') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        n = 0
        for row in readCSV:
            n += 1
    d = {
    'n':         n-1,
    'camera':    np.zeros((n-1,)),
    'frame':     np.zeros((n-1,)),
    'timestamp': np.zeros((n-1,)),
    'saturated': np.zeros((n-1,)),
    'rawmean':   np.zeros((n-1,)),
    'rawstd':    np.zeros((n-1,)),
    'mean':      np.zeros((n-1,)),
    'std':       np.zeros((n-1,)),
    'contrast':  np.zeros((n-1,)),
    'temperature':  np.zeros((n-1,))
    }

    # read data
    with open(fname) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        i = 0
        for row in readCSV:
            if i>0 and len(row)==10:
                d['camera'][i-1] = row[0]
                d['frame'][i-1] = row[1]
                d['timestamp'][i-1] = row[2]
                d['saturated'][i-1] = row[3]
                d['rawmean'][i-1] = row[4]
                d['rawstd'][i-1] = row[5]
                d['mean'][i-1] = row[6]
                d['std'][i-1] = row[7]
                d['contrast'][i-1] = row[8]
                d['temperature'][i-1] = row[9]
            elif i>0:
                print('CSV Format Has Changed!!!')
            i += 1

    return d

def sortbycamera_proto2(data, cids, nt):

    nc = len(cids)

    d = {
    'camera':    np.zeros((nt, nc)),
    'frame':     np.zeros((nt, nc)),
    'timestamp': np.zeros((nt, nc)),
    'saturated': np.zeros((nt, nc)),
    'rawmean':   np.zeros((nt, nc)),
    'rawstd':    np.zeros((nt, nc)),
    'mean':      np.zeros((nt, nc)),
    'std':       np.zeros((nt, nc)),
    'contrast':  np.zeros((nt, nc)),
    'temperature':  np.zeros((nt, nc))
    }

    ind = np.zeros((nc,))
    for i in range(data['n']):
        j = cids.index(data['camera'][i].astype('int').astype('str'))
        k = int(ind[j])
        d['camera'][k, j] = data['camera'][i]
        d['frame'][k, j] = data['frame'][i]
        d['timestamp'][k, j] = data['timestamp'][i]
        d['saturated'][k, j] = data['saturated'][i]
        d['rawmean'][k, j] = data['rawmean'][i]
        d['rawstd'][k, j] = data['rawstd'][i]
        d['mean'][k, j] = data['mean'][i]
        d['std'][k, j] = data['std'][i]
        d['contrast'][k, j] = data['contrast'][i]
        d['temperature'][k, j] = data['temperature'][i]
        ind[j] += 1

    return d

def removedark2(data, data_dict):

    mean = np.abs(data['rawmean'][1::2, :] - data['rawmean'][0::2, :])
    std = np.sqrt(np.abs(data['rawstd'][1::2, :]**2 - data['rawstd'][0::2, :]**2))
    mean_mean = np.mean(mean, axis=0)
    s = mean.shape
    for i in range(1,s[0]):
        for j in range(s[1]):
            if mean[i, j]<0.3*mean_mean[j]: ##********* clamping to a min value?
                mean[i, j] = mean[i-1, j]
                std[i, j] = std[i-1, j]
    contrast = std / mean

    # Calculates mean/std/contrast using dark frames from before and after (includes shot noise corrections)
    meanBA = np.abs(data['rawmean'][1::1, :] - data['rawmean'][0:-1:1, :])
    meanBA = (meanBA[1::2, :] + meanBA[0:-1:2, :])/2
    stdBA = np.sqrt(np.abs(data['rawstd'][1::1, :]**2 - data['rawstd'][0:-1:1, :]**2))
    stdBA = (stdBA[1::2, :] + stdBA[0:-1:2, :])/2
    # stdBA = np.sqrt(stdBA**2 - data_dict['K']*data_dict['gain']*meanBA)
    contrastBA = stdBA / meanBA

    d = {
    'camera':    data['camera'][1::2, :],
    'timestamp': data['timestamp'][1::2, :],
    'mean':      mean,
    'std':       std,
    'contrast':  contrast,
    'temperature':  data['temperature'][1::2, :],
    'mean_dark': data['rawmean'][0::2, :],
    'std_dark': data['rawstd'][0::2, :],
    'contrast_dark': data['contrast'][0::2, :], #fix
    'rawmean': data['rawmean'],
    'rawstd': data['rawstd'],
    'meanBA': meanBA,
    'stdBA': stdBA,
    'contrastBA': contrastBA
    }
    """
    d = {
    'camera':    data['camera'][1::2, :],
    'timestamp': data['timestamp'][1::2, :],
    'mean':      data['mean'][1::2, :],
    'std':       data['std'][1::2, :],
    'contrast':  data['contrast'][1::2, :],
    'temperature':  data['temperature'][1::2, :],
    'mean_dark': data['mean'][0::2, :],
    'std_dark': data['std'][0::2, :],
    'contrast_dark': data['contrast'][0::2, :],
    'rawmean': data['rawmean'],
    'rawstd': data['rawstd']
    }

    if np.sum(d['mean'])==0:

        d['mean'] = d['rawmean'][0::2, :] - d['rawmean'][1::2, :]
        var = d['rawstd'][0::2, :]**2 - d['rawstd'][1::2, :]**2
        var[var<0]=0
        d['std']=np.sqrt(var)
        d['contrast'] = d['std']/d['mean']
        d['mean_dark'] = d['rawmean'][1::2, :]
        d['std_dark'] = d['rawstd'][1::2, :]
        d['contrast_dark'] = d['std_dark']/d['mean_dark']
    """
    return d

def croptime_proto(d, data_dict):

    nstart = data_dict.get('nstart', 0)
    nstop = data_dict.get('nstop', d['mean'].shape[0]+1)
    for key in d.keys():
        d[key] = d[key][nstart:nstop, :]

    return d

def shotnoise(d, data_dict):

    d['std'] = np.sqrt( d['std']**2 - data_dict['K']*data_dict['gain']*d['mean'] )
    d['contrast'] = d['std']/d['mean']

    d['stdBA'] = np.sqrt( d['stdBA']**2 - data_dict['K']*data_dict['gain']*d['meanBA'] )
    d['contrastBA'] = d['stdBA']/d['meanBA']

    return d

def specialsauceInd(d):

    s = d['mean'].shape
    x = d['mean']#**2
    y = d['std']#**2

    sx = np.sum(x, axis=0)
    sy = np.sum(y, axis=0)
    sxx = np.sum(x**2, axis=0)
    sxy = np.sum(x*y, axis=0)
    a = (sxy - sx*sy/s[0])/(sxx - sx*sx/s[0])
    b = (sy - a*sx)/s[0]
    d['std'] -= b
    d['contrast'] = d['std']/d['mean']
    d['ncontrast_mean'] = a

    return d

def specialsauceAll(d):

    x = combinepositions(d, 'mean')
    y = combinepositions(d, 'std')
    n = x.shape[0]
    sx = np.sum(x, axis=0)
    sy = np.sum(y, axis=0)
    sxx = np.sum(x**2, axis=0)
    sxy = np.sum(x*y, axis=0)
    a = (sxy - sx*sy/n)/(sxx - sx*sx/n)
    b = (sy - a*sx)/n
    print('slope = ' + str(a) + ' offset = ' + str(b))
    for i in range(len(d)):
       d[i]['std'] -= b
       d[i]['contrast'] = d[i]['std']/d[i]['mean']

    return d

def removeOffsetAll(d, data_dict):

    x = combinepositions(d, 'mean')
    y = combinepositions(d, 'std')
    n = x.shape[0]
    nc = x.shape[1]

    if data_dict['offset']==0:
        a = np.mean(y/x, axis=0)
        b = np.zeros((nc,))
        c = np.zeros((nc,))
    elif data_dict['offset']==1:
        #1st way: normal equations for ax+b
        sx = np.sum(x, axis=0)
        sy = np.sum(y, axis=0)
        sxx = np.sum(x**2, axis=0)
        sxy = np.sum(x*y, axis=0)
        a = (sxy - sx*sy/n)/(sxx - sx*sx/n)
        b = (sy - a*sx)/n
        c = np.zeros((nc,))
        """
        print('slope = ' + str(a) + ' offset = ' + str(b))
        fig, ax = plt.subplots(nrows=nc, ncols=1)
        for i in range(nc):
            yf = np.zeros(x.shape)
            xf = np.array([np.amin(x[:, i]), np.amax(x[:, i])])
            yf = a[i]*xf+b[i]
            ax[i].plot(x[:, i], y[:, i], '.')
            ax[i].plot(xf, yf)
        d[0]['a_all'] = a
        d[0]['b_all'] = b
        d[0]['c_all'] = np.zeros((nc,))
        """
    elif data_dict['offset']==2:
        #2nd way: svd with 1/x term
        a = np.zeros((nc,))
        b = np.zeros((nc,))
        c = np.zeros((nc,))
        for i in range(nc):
            M = np.stack((x[:, i], np.ones(n,), 1/x[:, i]), axis=1)
            out = np.linalg.lstsq(M, y[:, i], rcond=None)
            if out[0][0]>1e10:
                a[i]=out[0][0]
                b[i]=out[0][1]
                c[i]=out[0][2]
            else:
                M = np.stack((x[:, i], np.ones(n,)), axis=1)
                out = np.linalg.lstsq(M, y[:, i], rcond=None)
                a[i]=out[0][0]
                b[i]=out[0][1]
                c[i]=0

    print('a = ' + str(a) + ' b = ' + str(b) + ' c = ' + str(c))
    col = data_dict['col']
    row = data_dict['row']
    cc = 'rgbcmkrgbcmkrgbcmkrgbcmkrgbcmk'
    rr = '<><><><><><>'
    yf = np.zeros(x.shape)
    x1 = np.amin(x, axis=0)
    x2 = np.amax(x, axis=0)

    fig, ax = plt.subplots(nrows=nc, ncols=1)
    for i in range(nc):
        for j in range(len(d)):
            ax[i].plot(d[j]['mean'][:, i], d[j]['std'][:, i], cc[col[j]]+rr[row[j]])
        xf = np.arange(np.floor(x1[i]), np.ceil(x2[i]))
        yf = a[i]*xf+b[i]+c[i]/xf
        ax[i].plot(xf, yf, color=[0.5, 0.5, 0.5])

    d[0]['a_all'] = a
    d[0]['b_all'] = b
    d[0]['c_all'] = c

    for i in range(len(d)):
        for j in range(nc):
            if data_dict['offset']==2:
                d[i]['std'][:, j] -= (b[j] + c[j]/d[i]['mean'][:, j])
            else:
                d[i]['std'][:, j] -= b[j]
        d[i]['contrast'] = d[i]['std']/d[i]['mean']

    return d

def specialsauceAlmostAll(d):


    x = combinepositions(d, 'mean')
    y = combinepositions(d, 'std')
    for i in range(x.shape[1]):
        m = np.stack((x[:, i], y[:, i]), axis=0)
        cov_mat = np.cov(m)
        eigen_vals, eigen_vecs = np.linalg.eig(cov_mat)

        print(cov_mat)
        print(eigen_vals)
        print(eigen_vecs)


    n = x.shape[0]
    sx = np.sum(x, axis=0)
    sy = np.sum(y, axis=0)
    sxx = np.sum(x**2, axis=0)
    sxy = np.sum(x*y, axis=0)
    a = (sxy - sx*sy/n)/(sxx - sx*sx/n)
    b = (sy - a*sx)/n
    print('slope = ' + str(a) + ' offset = ' + str(b))
    for i in range(len(d)):
       d[i]['std'] -= b
       d[i]['contrast'] = d[i]['std']/d[i]['mean']

    return d




def averagetimeseries_proto(d):

    d['mean_mean'] = np.nanmean(d['mean'], axis=0)
    d['std_mean'] = np.nanmean(d['std'], axis=0)
    d['contrast_mean'] = np.nanmean(d['contrast'], axis=0)
    d['mean_dark_mean'] = np.nanmean(d['mean_dark'], axis=0)
    d['std_dark_mean'] = np.nanmean(d['std_dark'], axis=0)
    d['contrast_dark_mean'] = np.nanmean(d['contrast_dark'], axis=0)

    return d

def getridofzeros(d):

    for key in d.keys():
        if np.array(d[key]).size == 4:
            for i in range(4):
                if np.array(d[key][i]).size==1:
                    if d[key][i]==0:
                        d[key][i]=np.nan

    return d

def combinepositions(d, keyname):

    npositions = len(d)  #number of positions
    if isinstance(d[0][keyname], bool) or isinstance(d[0][keyname], int) or isinstance(d[0][keyname], float):
        nt=1
        nc=1
    elif len(d[0][keyname].shape)==1:
        nt=1
        nc=d[0][keyname].shape[0]
    else:
        nt = d[0][keyname].shape[0] #number acquisition per position
        nc = d[0][keyname].shape[1] #number of cameras
    x = np.zeros((npositions*nt, nc))
    for i in range(npositions):
        m = int(i*nt)
        n = int((i+1)*nt)
        if nt==1:
            x[m:n, :] = np.expand_dims(d[i][keyname], axis=0)
        else:
            x[m:n, :] = d[i][keyname]

    return x

def combinepositions2(d, keyname, positions=0):

    npositions = len(d)  #number of positions
    if len(d[0][keyname].shape)==1:
        nt=1
        nc=d[0][keyname].shape[0]
    else:
        nt = d[0][keyname].shape[0] #number acquisition per position
        nc = d[0][keyname].shape[1] #number of cameras
    x = np.zeros((npositions*nt, nc))
    for i in range(npositions):
        m = int(i*nt)
        n = int((i+1)*nt)
        if nt==1:
            x[m:n, :] = np.expand_dims(d[i][keyname], axis=0)
        else:
            x[m:n, :] = d[i][keyname]

    return x


def plotscandataintime(d, data_dict):

    npositions = len(d)  #number of positions
    nt = d[0]['mean'].shape[0] #number acquisition per position
    nc = d[0]['mean'].shape[1] #number of cameras

    mean = combinepositions(d, 'mean')
    std =combinepositions(d, 'std')
    contrast = combinepositions(d, 'contrast')

    col = data_dict['col']
    row = data_dict['row']
    c = 'rgbcmkrgbcmkrgbcmkrgbcmkrgbcmk'
    r = '<><><><><>'
    """
    fig, ax = plt.subplots(nrows=nc, ncols=3)
    for j in range(npositions):
        for i in range(nc):
            ax[i, 0].plot(np.arange(j*nt, (j+1)*nt), d[j]['mean'][:, i], c[col[j]]+r[row[j]])
            ax[i, 1].plot(np.arange(j*nt, (j+1)*nt), d[j]['std'][:, i], c[col[j]]+r[row[j]])
            ax[i, 2].plot(np.arange(j*nt, (j+1)*nt), d[j]['contrast'][:, i], c[col[j]]+r[row[j]])
            ax[i, 0].plot(j*np.array([nt, nt]), np.array([np.amin(mean[:, i], axis=0), np.amax(mean[:, i], axis=0)]), ':r')
            ax[i, 1].plot(j*np.array([nt, nt]), np.array([np.amin(std[:, i], axis=0), np.amax(std[:, i], axis=0)]), ':r')
            ax[i, 2].plot(j*np.array([nt, nt]), np.array([np.amin(contrast[:, i], axis=0), np.amax(contrast[:, i], axis=0)]), ':r')
    fig.suptitle('Mean              Std              Contrast')
    """
    fig, ax = plt.subplots(nrows=nc, ncols=2)
    for j in range(npositions):
        for i in range(nc):
            ax[i, 0].plot(d[j]['mean'][:, i], d[j]['std'][:, i],      c[col[j]]+r[row[j]])
            ax[i, 1].plot(d[j]['mean'][:, i], d[j]['contrast'][:, i], c[col[j]]+r[row[j]])
            ax[i, 0].set_xlabel('Mean')
            ax[i, 0].set_ylabel('Std')
            ax[i, 1].set_xlabel('Mean')
            ax[i, 1].set_ylabel('Contrast')
    """
    fig, ax = plt.subplots(nrows=nc, ncols=1)
    for j in range(npositions):
        for i in range(nc):
            ax[i].plot(np.arange(j*nt, (j+1)*nt), d[j]['contrast'][:, i], 'k.-')
            ax[i].plot(j*np.array([nt, nt]), np.array([np.amin(contrast[:, i], axis=0), np.amax(contrast[:, i], axis=0)]), ':r')
    fig.suptitle('Contrast')
    """

def getoffsets(d):

    fig, ax = plt.subplots(nrows=1, ncols=4)
    for i in range(2):
        for j in range(4):
            ax[j].plot((d['std'][:, 0]-i)/d['mean'][:, j])


def getpulse(d, data_dict):

    nrho = len(data_dict['rho'])
    N = 512
    m = np.mean(d['contrast'], axis = 0)
    #f = np.fft.fftshift(np.fft.fft(d['contrast']-m, n=N, axis=0), axes=0)
    #d['acontrast'] = np.angle(f)
    d['pcontrast'] = np.abs(np.fft.fftshift(np.fft.fft(d['contrast']-m, n=N, axis=0), axes=0))
    p = d['pcontrast'][int(N/2):, :]
    d['pulse']=np.zeros((nrho,))
    d['amplitude']=np.zeros((nrho,))
    d['modulationdepth']=np.zeros((nrho,))
    ind = np.argmax(p, axis=0)
    for i in range(nrho):
        if ind[i]>24 and ind[i]<122:
            d['pulse'][i] = np.array([60*ind[i]/N/data_dict.get('dt')])
            #d['amplitude'][i] = p[ind[i], i]
            #d['modulationdepth'][i] = p[ind[i], i]/m[i]/100
        else:
            d['pulse'][i] = np.array([np.nan])
            #d['amplitude'][i] = np.array([np.nan])
            #d['modulationdepth'][i] = np.array([np.nan])

    return d

def highpass(d, data_dict):

    cut = 4
    #nrho = len(data_dict['rho'])
    #N = 512
    m = np.mean(d['contrast'], axis = 0)
    bob = np.fft.fft(d['contrast']-m, axis=0)
    bob[:cut, :] = 0
    bob[-cut:, :] = 0
    bob = np.real(np.fft.ifft(bob, axis=0))
    d['hcontrast'] = bob
    m = np.mean(d['mean'], axis = 0)
    bob = np.fft.fft(d['mean']-m, axis=0)
    bob[:cut, :] = 0
    bob[-cut:, :] = 0
    bob = np.real(np.fft.ifft(bob, axis=0))
    d['hmean'] = bob



    return d
"""
def getpulseplus(d, data_dict):

    nrho = len(data_dict['rho'])
    N = 512
    m = np.mean(d['contrast'], axis = 0)
    d['pcontrast'] = np.abs(np.fft.fftshift(np.fft.fft(d['contrast']-m, n=N, axis=0), axes=0))
    p = d['pcontrast'][int(N/2):, :]
    d['pulse']=np.zeros((nrho,))
    d['amplitude']=np.zeros((nrho,))
    d['modulationdepth']=np.zeros((nrho,))
    ind = np.argmax(p, axis=0)
    for i in range(nrho):
        if ind[i]>24 and ind[i]<122:
            d['pulse'][i] = np.array([60*ind[i]/N/data_dict.get('dt')])
            #d['amplitude'][i] = p[ind[i], i]
            #d['modulationdepth'][i] = p[ind[i], i]/m[i]/100
        else:
            d['pulse'][i] = np.array([np.nan])
            #d['amplitude'][i] = np.array([np.nan])
            #d['modulationdepth'][i] = np.array([np.nan])

    return d
"""

def getIntensityCalibration(params_dict):

    mu_eff = np.sqrt(3*params_dict['mua']*(params_dict['musp']+params_dict['mua']))
    G = G_semiinf(mu_eff, params_dict['rho'])
    ci = G/params_dict['I0'].T
    ci /= np.max(ci)

    return ci

def getContrastCalibration(params_dict):

    rhoplus = get_rhoplus_multi_c(params_dict)
    nexp = len(params_dict['exp'])
    exps = rhoplus[-nexp:]
    expinds = (np.round(exps/rhoplus[5])).astype('int32')
    rhoplus = rhoplus[:-nexp]
    calculated = specklecontrast_numerical_fitAlphaDb([params_dict['Db']], rhoplus)
    measured = params_dict['C0'] #- params_dict['dcontrast']
    cc = np.squeeze(calculated[expinds, :].T) / measured

    return cc

def G_semiinf(k, rho):

    G = np.exp(-k*rho)/(rho**2)

    return G

def get_rhoplus_multi_c(params_dict):

    nrho = len(params_dict['rho'])
    nexp = len(params_dict['exp'])

    rhoplus = np.zeros((8+nrho+nexp,))
    rhoplus[0] = params_dict['n'] # index of refraction
    rhoplus[1] = params_dict['wv'] # wavelength, mm
    rhoplus[2] = params_dict['mua'] # absorption coefficient, mm^-1
    rhoplus[3] = params_dict['musp'] # reduced scattering coefficient, mm^-1
    rhoplus[4] = params_dict['Tmax'] # maximum exposure time being simulated, s
    rhoplus[5] = params_dict['dtau'] # discretization of time step length, s
    rhoplus[6] = params_dict['zb'] # extrapolated boundary (for method of images), mm
    rhoplus[7] = params_dict['beta'] # constant from 0 to 1 depending on pixel size, polarization, and abberations, unitless
    rhoplus[8:(8+nrho)] = params_dict['rho'] # distance between source and detector fibers, mm
    rhoplus[(8+nrho):] = params_dict['exp']

    return rhoplus

def get_rhoplus_multi(params_dict):

    c2u = params_dict.get('cameras2use', np.arange(0, len(params_dict['rho'])+1))

    nrho = len(c2u)#len(params_dict['rho'])
    nexp = len(params_dict['exp'])

    rhoplus = np.zeros((8+nrho+nexp,))
    rhoplus[0] = params_dict['n'] # index of refraction
    rhoplus[1] = params_dict['wv'] # wavelength, mm
    rhoplus[2] = params_dict['mua'] # absorption coefficient, mm^-1
    rhoplus[3] = params_dict['musp'] # reduced scattering coefficient, mm^-1
    rhoplus[4] = params_dict['Tmax'] # maximum exposure time being simulated, s
    rhoplus[5] = params_dict['dtau'] # discretization of time step length, s
    rhoplus[6] = params_dict['zb'] # extrapolated boundary (for method of images), mm
    rhoplus[7] = params_dict['beta'] # constant from 0 to 1 depending on pixel size, polarization, and abberations, unitless
    rhoplus[8:(8+nrho)] = params_dict['rho'][c2u] # distance between source and detector fibers, mm
    rhoplus[(8+nrho):] = params_dict['exp']

    return rhoplus

def specklecontrast_numerical_fitAlphaDb(x, rhoplus):

    n = rhoplus[0]    # index of refraction
    wv = rhoplus[1]   # wavelength, mm
    mua = rhoplus[2]  # absorption coefficient, mm^-1
    musp = rhoplus[3] # reduced scattering coefficient, mm^-1
    Tmax = rhoplus[4] # maximum exposure time being simulated, s
    dtau = rhoplus[5] # discretization of time step length, s
    zb = rhoplus[6]   # extrapolated boundary (for method of images), mm
    alpha = 1         # fraction of scatterers that are moving
    if len(x)==1:
        beta = rhoplus[7]
        Db = x[0]     # diffusion coefficient of the scatterers, mm^2/s
        offset = 0
    elif len(x)==2:
        beta = x[0]   # constant from 0 to 1 depending on pixel size, polarization, and abberations, unitless
        Db = x[1]     # diffusion coefficient of the scatterers, mm^2/s
        offset = 0
    elif len(x)==3:
        beta = x[0]   # constant from 0 to 1 depending on pixel size, polarization, and abberations, unitless
        Db = x[1]     # diffusion coefficient of the scatterers, mm^2/s
        offset = x[2]
    rho = rhoplus[8:] # distance between source and detector fibers, mm

    ### calculations ###
    k0 = 2*np.pi*n/wv      #wavenumber
    tau = np.arange(0, Tmax+dtau, dtau)
    dr2 = (0.5*(tau[1:]+tau[:-1]))*(6*Db)
    z0 = 1/musp
    contrast = np.zeros((len(dr2), len(rho)))   # (exposures, rhos)
    r1 = np.sqrt( rho**2 + z0**2 )                 # len = rhos
    r2 = np.sqrt( rho**2 + (z0+(2*zb))**2 )        # len = rhos
    k = np.sqrt( alpha*(k0**2)*(musp**2)*dr2 + 3*mua*musp ) # len = exposures
    sourceterm = np.exp(-np.outer(k,r1))/np.tile(r1,(len(dr2), 1))
    imageterm =   np.exp(-np.outer(k,r2))/np.tile(r2,(len(dr2), 1))
    G1 = sourceterm - imageterm
    g1 = np.zeros(G1.shape)  # (exposures, rhos)
    for i in range(len(rho)):
        g1[:, i] = G1[:, i]/G1[0,i]
    taus = np.tile(tau[:,np.newaxis],(1, len(rho)))
    part1 = np.cumsum((g1**2), axis=0) / taus[1:]
    part2 = np.cumsum(0.5*(taus[1:]+taus[:-1])*(g1**2), axis=0) / (taus[1:]**2)
    part1minuspart2 = part1 - part2
    contrast = np.sqrt( 2 * beta * dtau * (part1minuspart2) ) + offset

    return contrast

def fitattenuation(mean, p, pp, musp, sigma, displayfits):

    nt = mean.shape[0]
    nexp = mean.shape[2]

    ### get initial guess
    logdata = np.log((p**2)*mean[0, :, 0])
    M = M = np.stack((np.ones(len(p),), -p), axis=1)
    Minv = np.linalg.pinv(M)
    temp = Minv @ logdata
    initial_guess = np.array([np.exp(temp[0]), temp[1]])

    I0 = np.zeros((nt, nexp))
    if len(initial_guess)==2:
        fitI0 = True
    elif len(initial_guess)==1:
        fitI0 = False

    mu_eff = np.zeros((nt, nexp))
    mua = np.zeros((nt, nexp))

    mean_calc = np.zeros((nt, len(pp), nexp))
    for i in range(nt):
        for j in range(nexp):
            data = np.squeeze(mean[i, :, j])
            uncertainty = np.squeeze(sigma[i, :, j])
            if np.all(np.isfinite(data)) and np.all(np.isfinite(uncertainty)):
                out = least_squares(residual_intensity3, initial_guess, bounds=(0, np.inf), args=(p, data, uncertainty))
                if fitI0==1:
                    mu_eff[i, j] = out.x[1]
                    I0[i, j] = out.x[0]
                else:
                    mu_eff[i, j] = out.x[0]
                    I0[i, j] = 1
                mua[i, j] = mu_eff[i, j]**2/(3*musp[j])
                mean_calc[i, :, j] = I0[i, j]*G_semiinf(mu_eff[i, j], pp)
            else:
                mu_eff[i, j] = np.nan
                I0[i, j] = np.nan
                mua[i, j] = np.nan
                mean_calc[i, :, j] = np.nan

    if displayfits:
        colors = 'rbgcmykrgbcmykrgbcmykrgbcmykrbgcmykrgbcmykrgbcmykrgbcmykrbgcmykrgbcmykrgbcmykrgbcmyk'
        """
        num1 = 0
        num2 = 0
        plt.figure()
        for i in range(0, nt):
            for j in range(nexp):
                plt.plot(pp, mean_calc[i, :, j], color=colors[num1])
                num1 += 1
            for j in range(nexp):
                plt.plot(p, np.squeeze(mean[i, :, j]), '*', color=colors[num2])
                num2 += 1
        plt.xlabel('Distance (mm)')
        plt.ylabel('Intensity (ADU)')
        plt.title('u_eff = ' + np.array2string(np.squeeze(mu_eff), precision = 2) + ' mm-1')
        plt.legend(('fit', 'data'))
        """
        num1 = 0
        num2 = 0
        plt.figure()
        for i in range(0, nt):
            for j in range(nexp):
                plt.semilogy(pp, (pp**2)*mean_calc[i, :, j], color=colors[num1])
                num1 += 1
            for j in range(nexp):
                plt.semilogy(p, (p**2)*np.squeeze(mean[i, :, j]), '*', color=colors[num2])
                num2 += 1
        plt.xlabel('Distance (mm)')
        plt.ylabel('Distance^2 * Intensity (mm^2 * ADU)')
        plt.title('u_eff')# + np.array2string(np.squeeze(mu_eff), precision = 2) + ' mm-1')
        plt.legend(('fit', 'data'))

    return mu_eff, mua

def residual_intensity3(x, rho, data, sigma):

    if len(x)==2:
        R = (x[0]*G_semiinf(x[1], rho) - data)/sigma
    elif len(x)==1:
        R = (G_semiinf(x[0], rho) - data)/sigma

    return R

def residual_intensity2(x, rho, data):

    if len(x)==2:
        R = x[0]*G_semiinf(x[1], rho) - data
    elif len(x)==1:
        R = G_semiinf(x[0], rho) - data
        #print(R)
    return R

def fitflow_multi(params_dict, contrast, p, pp, exp, expp, mua):

    nt = contrast.shape[0]

    fitbeta = params_dict.get('fitbeta', 0)
    if fitbeta==2:
        initial_guess = [0.3, 1e-6, 0]    #[beta, alpha*Db, offset]
        lb = [0, 0, 0]
        ub = [1, np.inf, np.inf]
    elif fitbeta==1:
        initial_guess = [0.3, 1e-6]    #[beta, alpha*Db]
        lb = [0, 0]
        ub = [1, np.inf]
    else:
        initial_guess = [1e-6] #[alphaDb]
        lb = 0
        ub = np.inf

    flowseparate = params_dict.get('flowseparate', False)
    if flowseparate==True:
        alphaDb = np.zeros((nt, len(p)))
        beta = np.zeros((nt, len(p)))
        offset = np.zeros((nt, len(p)))
        new_dict = copy.deepcopy(params_dict)
    else:
        alphaDb = np.zeros((nt, 1))
        beta =  np.zeros((nt, 1))
        offset = np.zeros((nt, 1))
        rhoplus = get_rhoplus_multi(params_dict)

    contrast_calc = np.zeros((nt, len(pp), len(expp)))

    for i in range(nt):
        for j in range(alphaDb.shape[1]):
            if params_dict.get('usemueff', 1)==1:
                params_dict['mua']=np.mean(mua[i])
                rhoplus = get_rhoplus_multi(params_dict)
            if flowseparate == True:
                data = np.expand_dims(contrast[i, j, :], axis=0)
                new_dict['rho'] = [params_dict['rho'][j]]
                rhoplus = get_rhoplus_multi(new_dict)
            else:
                data = contrast[i, :, :]

            if np.all(np.isfinite(data)) and np.all(np.isfinite(rhoplus)):
                out = least_squares(residual_speckle_multi, initial_guess, bounds=(lb , ub), args=(rhoplus, data))

                if fitbeta==2:
                    beta[i, j] = out.x[0]
                    alphaDb[i, j] = out.x[1]
                    offset[i, j] = out.x[2]
                if fitbeta==1:
                    beta[i, j] = out.x[0]
                    alphaDb[i, j] = out.x[1]
                else:
                    beta[i, j] = params_dict['beta']
                    alphaDb[i, j] = out.x[0]
            else:
                beta[i, j] = np.nan
                alphaDb[i, j] = np.nan

            if flowseparate==True:
                rhoplus2 = np.append(rhoplus[:8], np.append(p[j], expp))
                expinds = (np.round(expp/rhoplus[5])).astype('int32')
                temp = specklecontrast_numerical_fitAlphaDb(out.x, rhoplus2[:-len(expp)])
                contrast_calc[i, j, :] = temp[expinds, :].T
            else:
                rhoplus2 = np.append(rhoplus[:8], np.append(pp, expp))
                expinds = (np.round(expp/rhoplus[5])).astype('int32')
                temp = specklecontrast_numerical_fitAlphaDb(out.x, rhoplus2[:-len(expp)])
                contrast_calc[i, :, :] = temp[expinds, :].T

    colors = 'rbgcmykrgbcmykrgbcmykrgbcmykrbgcmykrgbcmykrgbcmykrgbcmykrbgcmykrgbcmykrgbcmykrgbcmyk'
    num1 = 0
    num2 = 0
    nexp = len(params_dict['exp'])
    if params_dict['displayfits'] and nt<40:
        plt.figure()
        for i in range(0, nt):
            for j in range(nexp):
                plt.plot(pp, contrast_calc[i, :, j], color=colors[num1])
                num1 += 1
            for j in range(nexp):
                plt.plot(p, np.squeeze(contrast[i, :, j]), '*', color=colors[num2])
                num2 += 1
        plt.xlabel('Distance (mm)')
        plt.ylabel('Contrast')
        plt.title('Blood Flow Index')# = ' + np.array2string(alphaDb, precision = 2) + ' mm2/s')
        plt.legend(('fit', 'data'))

    return alphaDb, beta, offset

def residual_speckle_multi(x, rhoplus, data):

    exps = rhoplus[-data.shape[1]:]
    expinds = (np.round(exps/rhoplus[5])).astype('int32')
    rhoplus = rhoplus[:-data.shape[1]]
    calculated = specklecontrast_numerical_fitAlphaDb(x, rhoplus)
    R = calculated[expinds, :] - data.T

    return R.flatten()

def plotlocation(data, data_dict, title):

    nt = data[title].shape[0]
    nc = data[title].shape[1]
    dt = data_dict['dt']
    t = np.arange(0, dt*nt, dt)
    if title=='pcontrast':
        t=np.arange(nt)
    fig, ax = plt.subplots(nrows=nc, ncols=1)
    ax[0].set_title(title)
    ax[nc-1].set_xlabel('seconds')
    for j in range(nc):
        if np.isfinite(data['pulse'][j]) == False:
            c=[1.0, 0.5, 0.25]
        else:
            c=[0, 0, 0]
        ax[j].plot(t, data[title][:, j], '.-', color=c)
        ax[j].set_ylabel(data_dict.get('cids')[j])

def plotlocationUpsidedown(data, data_dict, title):

    nt = data[title].shape[0]
    nc = data[title].shape[1]
    dt = data_dict['dt']
    t = np.arange(0, dt*nt, dt)
    fig, ax = plt.subplots(nrows=nc, ncols=1)
    ax[0].set_title('Upside down ' + title)
    ax[nc-1].set_xlabel('seconds')
    for j in range(nc):
        p = copy.deepcopy(data[title][:, j])
        p -= np.amin(p, axis=0)
        p = 1 - p
        #if np.isfinite(data['pulse'][j]) == False:
        #    c=[1.0, 0.5, 0.25]
        #else:
        c=[0, 0, 0]
        ax[j].plot(t, p, '.-', color=c)
        ax[j].set_ylabel(data_dict.get('cids')[j])
    """
    fig, ax = plt.subplots(nrows=nc, ncols=1)
    ax[0].set_title('Upside down ' + title)
    ax[nc-1].set_xlabel('seconds')
    for j in range(nc):
        p = data[title][:, j]
        #p /= np.amax(p, axis=0)
        p = 1/p
        if np.isfinite(data['pulse'][j]) == False:
            c=[1.0, 0.5, 0.25]
        else:
            c=[0, 0, 0]
        ax[j].plot(t, p, '.-', color=c)
        ax[j].set_ylabel(data_dict.get('cids')[j])
    """

def plotscandata(data, data_dict, title):

    nc = data[0][title].shape[1]
    row = data_dict.get('row')
    col = data_dict.get('col')
    ns = len(data)

    maxi = np.sort(data[0][title], axis=0)[-2, :]
    mini = np.sort(data[0][title], axis=0)[1, :]
    for i in range(1, ns):
        maxi = np.maximum(maxi, np.sort(data[i][title], axis=0)[-2, :])
        mini = np.minimum(mini, np.sort(data[i][title], axis=0)[1, :])

    fig, ax = plt.subplots(nrows=int(2*nc+1), ncols=int(ns/2))
    for i in range(ns):
        for j in range(nc):
            if data[i]['completedscan'] == False:
                c=[1, 0, 0]
            elif np.isfinite(data[i]['pulse'][j]) == False:
                c=[1.0, 0.5, 0.25]
            else:
                c=[0, 0, 0]
            ax[int(row[i]*int(nc+1)+j), int(col[i])].plot(data[i][title][:, j], '.-', color=c)
            if title == 'mua':
                ax[int(row[i]*int(nc+1)+j), int(col[i])].set_ylim([mini[j], maxi[j]])
    fig.suptitle(title)
    for i in range(int(ns/2)):
        ax[nc, i].axis('off')


def plotscan_slider(data, data_dict, selection):

    data = copy.deepcopy(data)
    nc = data[0][selection].shape[1]
    ns = len(data)

    def update_selection(change):
        global global_maxs
        global global_mins

        data_concat = np.concatenate([loc[data_selector.value] for loc in data], axis=0)

        # find min and max for each concatenated column (each camera)
        global_maxs = [max(data_concat[:,cam]) for cam in range(nc)]
        global_mins = [min(data_concat[:,cam]) for cam in range(nc)]
        update_values(None)


    def update_values(change):

        loc = loc_slider.value - 1

        fig.suptitle(f'{data_selector.value} - location {loc+1}')

        for cam in range(nc):

            if data[loc]['completedscan'] == False:
                c = [1, 0, 0]
            elif np.isfinite(data[loc]['pulse'][cam]) == False:
                c = [1.0, 0.5, 0.25]
            else:
                c = [0, 0, 0]

            lines[cam].set_ydata(data[loc][data_selector.value][:,cam])
            lines[cam].set_color(c)
            if scale_checkbox.value == True:
                axs[cam].set_ylim(global_mins[cam], global_maxs[cam])
                axs[cam].relim()
            else:
                axs[cam].autoscale()
                axs[cam].relim()
        fig.canvas.draw()


    loc_slider = widgets.IntSlider(
        value=1, min=1, max=30, orientation='horizontal', continuous_update=False, description='location')
    loc_slider.observe(update_values, 'value')

    scale_checkbox = widgets.Checkbox(value=False, description='global scaling')
    scale_checkbox.observe(update_values, 'value')

    data_options = ['camera', 'timestamp', 'mean', 'std', 'contrast', 'temperature', 'mean_dark', 'std_dark', 'contrast_dark', 'rawmean', 'rawstd', 'hcontrast', 'hmean']
    data_selector =  widgets.Dropdown(options=data_options, value=selection, description='value to plot')
    data_selector.observe(update_selection, 'value')

    output = widgets.Output()
    with output:

        fig, axs = plt.subplots(nc,1, figsize=(8,8))
        lines = [None] * nc
        for cam in range(nc):
            lines[cam], = axs[cam].plot(data[0]['std'][:,cam])

    update_selection(None)

    box = widgets.VBox( layout={'border': '1px solid black',} )
    box.children = [output, loc_slider, scale_checkbox, data_selector]
    display(box)
    return box


def plottogether(data, data_dict, title, sign):

    normalize = True
    plotStds = True
    nc = len(data[0][title])
    row = data_dict.get('row')
    col = data_dict.get('col')
    dt = data_dict.get('dt')
    ns = len(data)

    if data_dict.get('scanorder')=='hhc30':
        col = [0, 1, 2, 4, 5, 6,  2, 3, 4,  0, 1, 2, 4, 5, 6,  0, 1, 2, 4, 5, 6,  2, 3, 4,  0, 1, 2, 4, 5, 6]
        row = [0, 0, 0, 0, 0, 0,  3, 3, 3,  6, 6, 6, 6, 6, 6,  1, 1, 1, 1, 1, 1,  4, 4, 4,  7, 7, 7, 7, 7, 7]

    maxi = 0
    mini = 0
    for i in range(ns):
        for j in range(nc):
            if (i==0 and j==0) and all(np.isfinite(data[0][title][0])):
                maxi = np.amax(sign*data[0][title][0])
                mini = np.amin(sign*data[0][title][0])
            elif all(np.isfinite(data[i][title][j])):
                maxi = np.maximum(maxi, np.amax(-data[i][title][j]))
                mini = np.minimum(mini, np.amin(-data[i][title][j]))
    fig, ax = plt.subplots(nrows=np.amax(row)+1, ncols=np.amax(col)+1)
    c = 'rgbkcm'
    for i in range(ns):
        for j in range(nc):
            t = np.arange(0, len(data[i][title][j])*dt, dt)
            if normalize:
                f = sign*data[i][title][j]
                f_std = sign*data[i][title + '_std'][j]
                if all(np.isfinite(f)):
                    mini = np.amin(f)
                    f -= mini
                    maxi = np.amax(f)
                    f /= maxi
                    f_std /= maxi
                    ax[int(row[i]), int(col[i])].plot(t, f, '.-', color=c[j])
                    if plotStds:
                        ax[int(row[i]), int(col[i])].plot(t, f - f_std, linewidth=1, alpha=0.3, color=c[j])
                        ax[int(row[i]), int(col[i])].plot(t, f + f_std, linewidth=1, alpha=0.3, color=c[j])

            else:
                ax[int(row[i]), int(col[i])].plot(t, sign*data[i][title][j], '.-', color=c[j])
                if plotStds:
                    ax[int(row[i]), int(col[i])].plot(t, f - f_std, linewidth=1, alpha=0.3, color=c[j])
                    ax[int(row[i]), int(col[i])].plot(t, f + f_std, linewidth=1, alpha=0.3, color=c[j])
                ax[int(row[i]), int(col[i])].set_ylim((mini, maxi))

    if data_dict.get('scanorder')=='hhc30':
        for i in range(7):
            ax[2, i].set_axis_off()
            ax[5, i].set_axis_off()
        for i in range(8):
            if i!=3 and i!=4:
                ax[i, 3].set_axis_off()
            if i==3 or i==4:
                ax[i, 0].set_axis_off()
                ax[i, 1].set_axis_off()
                ax[i, 5].set_axis_off()
                ax[i, 6].set_axis_off()
        ax[0, 0].set_ylabel('Left')
        ax[1, 0].set_ylabel('Right')
        ax[3, 2].set_ylabel('Left')
        ax[4, 2].set_ylabel('Right')
        ax[6, 0].set_ylabel('Left')
        ax[7, 0].set_ylabel('Right')
        ax[0, 1].set_title('Forehead Verticle')
        ax[0, 5].set_title('Forehead Horizontal')
        ax[3, 3].set_title('Sylvian Fissure Diagonal')
        ax[6, 1].set_title('Temple Verticle')
        ax[6, 5].set_title('Temple Horizontal')
        for i in range(7):
            ax[5, i].set_xlabel('Seconds')
            for j in range(5):
                #ax[j, i].set_xticks([])
                if i>0:
                    ax[j, i].set_yticks([])

    fig.suptitle(title)
    if plotStds:
        fig.suptitle(title + ' with standard deviation')




def plot2scandata(data, data_dict, title1, title2, normalize, offset):

    row = data_dict.get('row')
    col = data_dict.get('col')
    ns = len(data)
    fig, ax = plt.subplots(nrows=9, ncols=int(ns/2))
    for i in range(ns):
        d1 = np.copy(data[i][title1])
        d2 = np.copy(data[i][title2])
        if offset==True:
            d1 -= np.amin(d1, axis=0)
            d2 -= np.amin(d2, axis=0)
            d1 /= np.amax(d1, axis=0)
            d2 /= np.amax(d2, axis=0)
        if normalize==True:
            d1 /= np.mean(data[i][title1], axis=0)
            d2 /= np.mean(data[i][title2], axis=0)

        ax[int(row[i]*5+0), int(col[i])].plot(d1[:, 0], '.-b')
        ax[int(row[i]*5+1), int(col[i])].plot(d1[:, 1], '.-b')
        ax[int(row[i]*5+2), int(col[i])].plot(d1[:, 2], '.-b')
        ax[int(row[i]*5+3), int(col[i])].plot(d1[:, 3], '.-b')
        ax[int(row[i]*5+0), int(col[i])].plot(d2[:, 0], '.-g')
        ax[int(row[i]*5+1), int(col[i])].plot(d2[:, 1], '.-g')
        ax[int(row[i]*5+2), int(col[i])].plot(d2[:, 2], '.-g')
        ax[int(row[i]*5+3), int(col[i])].plot(d2[:, 3], '.-g')
    fig.suptitle('blue = ' + title1 + ',   green = ' + title2)


def getdifference(data, title, col, row):

    if isinstance(data[0][title], bool) or isinstance(data[0][title], int) or isinstance(data[0][title], float):
        s=1
    else:
        s = data[0][title].shape[0]
    n = len(col)
    d = np.zeros((2, int(n/2), s))
    for i in range(n):
        d[row[i], col[i], :] = data[i][title]

    sub = d[0, :, :] - d[1, :, :]
    add = d[0, :, :] + d[1, :, :]
    normdiff = 2*sub/add

    return normdiff, sub, d

def showAverageData(data, data_dict, title):

    theta = data_dict.get('theta')
    Nx = 200
    Ny = 100
    h=len(data)
    hh=int(h/2)
    ns = len(data)
    a = Nx/2.2
    b = Ny/3
    x = a*np.cos(theta)+ Nx/2
    y = b*np.sin(theta)+ Ny/2

    normdiff, diff, d = getdifference(data, title, data_dict['col'],  data_dict['row'])

    if data_dict.get('scanorder')=='hhc30':
        fig, ax = plt.subplots(nrows=2, ncols=1)
        c = 'rgbk'
        ax[0].set_title(title)
        for i in range(d.shape[2]):
            for j in range(5):
                ax[0].plot([3*j, 3*j+1, 3*j+2], d[0, 3*j:3*(j+1), i], '<-', color=c[i])
                ax[0].plot([3*j, 3*j+1, 3*j+2], d[1, 3*j:3*(j+1), i], '>-', color=c[i])
        ax[0].set_ylabel(title)
        ax[0].set_xticks([])

        for j in range(5):
            for i in range(diff.shape[1]):
                ax[1].plot([3*j, 3*j+1, 3*j+2], diff[3*j:3*(j+1), i], '.-', color=c[i])
        ax[1].legend(data_dict['rho'])
        ax[1].set_title('Left Right Difference')
        ax[1].set_ylabel('Left - Right')#('2 * (L-R)/(L+R)')
        ax[1].set_xticks([1, 4, 7, 10, 13])
        ax[1].set_xticklabels(['Forehead V', 'Forehead H', 'Fissure D', 'Temple V', 'Temple H'])
        return
    elif data_dict.get('scanorder')=='generic30':
        d = np.concatenate((d[0, :, :], d[1, :, :]))
        groupNameList = data_dict['groupNameList']
        xAxisLabelLocations = list(np.arange((30/len(groupNameList)-1)/2,30,30/len(groupNameList)))
        smplPerGrp = int(30/len(groupNameList)) # samples per group, set to specific value if needed
        plotIndexs = np.array(list(range(smplPerGrp)))
        fig, ax = plt.subplots(nrows=2, ncols=1)
        c = 'rgbk'
        ax[0].set_title(title)
        for j in range(len(groupNameList)):
            for i in range(d.shape[1]):
                ax[0].plot(plotIndexs+smplPerGrp*j, d[smplPerGrp*j:smplPerGrp*(j+1), i], '<-', color=c[i])
        ax[0].set_ylabel(title)
        if 'legendNameList' in data_dict.keys():
            ax[0].legend(data_dict['legendNameList'])
        else:
            ax[0].legend(data_dict['rho'])
        ax[0].set_xticks(xAxisLabelLocations)
        ax[0].set_xticklabels(groupNameList)

        for j in range(len(groupNameList)):
            ax[1].plot(plotIndexs+smplPerGrp*j, d[smplPerGrp*j:smplPerGrp*(j+1), 1] - d[smplPerGrp*j:smplPerGrp*(j+1), 0], '<-', color=c[1])
        ax[1].set_ylabel('Cam2 minus Cam1')
        ax[1].set_xticks(xAxisLabelLocations)
        ax[1].set_xticklabels(groupNameList)
        return

    fig, ax = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios': [2, 1, 1]})
    fig.suptitle(title)
    im = np.nan*np.zeros((Ny, Nx))
    for i in range(ns):
        xx = int(x[i])
        yy = int(y[i])

        im[yy-h:yy-hh, xx-h:xx+h] = data[i][title][0]
        im[yy-hh:yy, xx-h:xx+h] = data[i][title][1]
        im[yy:yy+hh, xx-h:xx+h] = data[i][title][2]
        im[yy+hh:yy+h, xx-h:xx+h] = data[i][title][3]

        ax[0].text(xx-h, yy-h+3, np.array2string(data[i][title][0], precision=2), color=[1, 1, 1])
        ax[0].text(xx-h, yy-hh+3, np.array2string(data[i][title][1], precision=2), color=[1, 1, 1])
        ax[0].text(xx-h, yy+3, np.array2string(data[i][title][2], precision=2), color=[1, 1, 1])
        ax[0].text(xx-h, yy+hh+3, np.array2string(data[i][title][3], precision=2), color=[1, 1, 1])
    ax[0].imshow(im)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_xlabel('<- Back            Front ->')
    ax[0].set_ylabel('<- Right            Left ->')

    c = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']#'rgbcmk'
    for i in range(d.shape[2]):
        ax[1].plot(d[0, :, i], '<-', color=c[i])
        ax[1].plot(d[1, :, i], '>-', color=c[i])
    ax[1].set_ylabel(title)
    ax[1].set_xticks([])

    ax[2].plot(diff, '.-')
    ax[2].legend(data_dict['rho'])
    ax[2].set_title('Left Right Difference')
    ax[2].set_xlabel('<-- Back            Front -->')
    ax[2].set_ylabel('Left - Right')#('2 * (L-R)/(L+R)')
    ax[2].set_xticks([])

def showResult(data, data_dict, title):

    theta = data_dict.get('theta')
    Nx = 200
    Ny = 100
    h=len(data)
    hh=int(h/2)
    ns = len(data)
    a = Nx/2.2
    b = Ny/3
    x = a*np.cos(theta)+ Nx/2
    y = b*np.sin(theta)+ Ny/2

    normdiff, diff, d = getdifference(data, title, data_dict['col'],  data_dict['row'])

    if data_dict.get('scanorder')=='hhc30':
        fig, ax = plt.subplots(nrows=2, ncols=1)
        ax[0].set_title(title)
        for j in range(5):
            ax[0].plot([3*j, 3*j+1, 3*j+2], d[0, 3*j:3*(j+1), :], '<-k')
            ax[0].plot([3*j, 3*j+1, 3*j+2], d[1, 3*j:3*(j+1), :], '>-k')
        ax[0].set_ylabel(title)
        ax[0].set_xticks([])

        for j in range(5):
            ax[1].plot([3*j, 3*j+1, 3*j+2], diff[3*j:3*(j+1), :], '.-k')
        ax[1].legend(data_dict['rho'])
        ax[1].set_title('Left Right Difference')
        ax[1].set_ylabel('Left - Right')#('2 * (L-R)/(L+R)')
        ax[1].set_xticks([1, 4, 7, 10, 13])
        ax[1].set_xticklabels(['Forehead V', 'Forehead H', 'Fissure D', 'Temple V', 'Temple H'])
        return
    elif data_dict.get('scanorder')=='generic30':
        d = np.concatenate((d[0, :, :], d[1, :, :]))
        groupNameList = data_dict['groupNameList']
        xAxisLabelLocations = list(np.arange((30/len(groupNameList)-1)/2,30,30/len(groupNameList)))
        smplPerGrp = int(30/len(groupNameList)) # samples per group, set to specific value if needed
        plotIndexs = np.array(list(range(smplPerGrp)))
        fig, ax = plt.subplots(nrows=2, ncols=1)
        c = 'rgbk'
        ax[0].set_title(title)
        for j in range(len(groupNameList)):
            ax[0].plot(plotIndexs+smplPerGrp*j, d[smplPerGrp*j:smplPerGrp*(j+1), :], '<-', color='k')
        ax[0].set_ylabel(title)
        if 'legendNameList' in data_dict.keys():
            ax[0].legend(data_dict['legendNameList'])
        else:
            ax[0].legend(data_dict['rho'])
        ax[0].set_xticks(xAxisLabelLocations)
        ax[0].set_xticklabels(groupNameList)
        return

    fig, ax = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios': [2, 1, 1]})
    fig.suptitle(title)
    im = np.nan*np.zeros((Ny, Nx))
    for i in range(ns):
        xx = int(x[i])
        yy = int(y[i])
        im[yy-h:yy+h, xx-h:xx+h] = data[i][title][0]
        if np.isfinite(data[i][title]):
            value = np.array2string(data[i][title][0], precision=2)
        else:
            value = 'nan'
        ax[0].text(xx-hh*1.5, yy, value, color=[1, 1, 1])

    ax[0].imshow(im)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_xlabel('<- Back            Front ->')
    ax[0].set_ylabel('<- Right            Left ->')

    ax[1].plot(d[0, :, :], '<-k')
    ax[1].plot(d[1, :, :], '>-k')
    ax[1].set_ylabel(title)
    ax[1].set_xticks([])

    ax[2].plot(normdiff, '.-k')
    ax[2].set_title('Left Right Difference')
    ax[2].set_xlabel('<-- Back            Front -->')
    ax[2].set_ylabel('2 * (L-R)/(L+R)')
    ax[2].set_xticks([])


def getdatadictionary(scanorder, calname, scanname, scannername):

    if scanorder == 'zigzag':
        col = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5])
        row = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1])
        theta = 3.14*np.array([-0.85, 0.85, -0.7, 0.7, -0.55, 0.55, -0.42, 0.42, -0.28, 0.28, -0.12, 0.12])
    elif scanorder == 'circle':
        col = np.array([0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0])
        row = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1])
        theta = 3.14*np.array([-0.85, -0.7, -0.55, -0.42, -0.28, -0.12, 0.12, 0.28, 0.42, 0.55, 0.7, 0.85])
    elif scanorder == 'hhc30':
        #col = np.array([0, 1, 2, 0, 1, 2, 3, 4, 5, 3, 4, 5, 6, 7, 8, 6, 7, 8, 9, 10, 11, 9, 10, 11, 12, 13, 14, 12, 13, 14])
        #row = np.array([0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1])
        col = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        row = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        theta = np.nan
    elif scanorder == 'generic30':
        col = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        row = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        theta = np.nan
    data_dict = {
        'calname':  calname,
        'dataname': scanname,
        'col': col,
        'row': row,
        'theta': theta,
        'npositions': len(col),
        'scannername': scannername,
        #'nstart': 0,
        #'nstop': 40,
        'K': .12
        }

    return data_dict

def getparamsdictionary(batch):

    if batch:
        displayfits = False
    else:
        displayfits = True

    params_dict = {
    'n': 1.4,      # index of refraction
    'wv': np.array([8.5e-4]),    # wavelength, mm
    'mua': np.array([0.09]),   # absorption coefficient, mm^-1
    'musp': np.array([1.2]),     # reduced scattering coefficient, mm^-1
    'Db': 1e-6,    # diffusion coefficient of the scatterers, mm^2/s
    'dV2': 1,      # mean squared velocity of the scatterers, (mm/s)^2
    'Tmax': 1e-3,  # maximum exposure time being simulated, s
    'dtau': 1e-7,  # discretization of time step length, s
    'zb': 0.1,     # extrapolated boundary
    'cameras2use': np.array([0, 1, 2, 3]),
    'fitbeta': 0,
    'beta': 0.25,
    'displayfits': displayfits
    #'dcontrast': np.array([0.4, 0, 0, 0]),
    #'betas': np.array([0.25, 0.25, 0.25, 0.5])
    }

    return params_dict

def getjsoninfo(data_dict, params_dict):

    ### load json file ###
    with open(data_dict.get('dataname') + '/scan_metadata.json') as json_file:
        json_dict = json.load(json_file)
    nrho   = json_dict.get('cameraParameters').get('numCameras')
    scannerCameras = list(json_dict.get('cameraParameters').get('cameraInfo').keys())
    for i in range(len(scannerCameras)):
        if scannerCameras[i] == data_dict.get('scannername'):
            cids = list(json_dict.get('cameraParameters').get('cameraInfo').keys())[i+1:i+5]
            data_dict['cids'] = cids
    gain = np.zeros((nrho,))
    rho = np.zeros((nrho,))
    for i in range(nrho):
        gain[i]=json_dict.get('cameraParameters').get('cameraInfo').get(cids[i]).get('gain')
        rho[i]=json_dict.get('cameraParameters').get('cameraInfo').get(cids[i]).get('separation_mm')
    data_dict['gain'] = gain
    data_dict['rho'] = rho
    params_dict['rho'] = rho
    data_dict['nt']   = json_dict.get('cameraParameters').get('numImages')
    exp = np.array([json_dict.get('delayParameters').get('pulseWidth_s')])
    params_dict['exp']=exp
    #nexp = len(exp)
    data_dict['dt'] = 2/json_dict.get('cameraParameters').get('frameAcquisitionRate_Hz')
    data_dict['notes'] = json_dict.get('sampleParameters').get('experimentNotes')

    return data_dict, params_dict

def goldenpulse(data, data_dict, title):

    ### constants (that should be in data_dict) ###
    N=512        # zero padded length for FFT
    minbpm = 45  # minimum heartbeats per minute
    maxbpm = 120 # maximum heartbeats per minute

    ### data ###
    y = data[title]
    nc = y.shape[1]
    dt = data_dict['dt']
    data[title + '_waveform'] = nc*[0]
    data[title + '_waveform_std'] = nc*[0]
    data[title + '_waveform_std_mean'] = np.array(np.zeros((nc,)))
    data[title + '_waveform_passedCount'] = np.array(np.zeros((nc,)))

    ### approximate period ###
    m=np.mean(y, axis=0)
    Y=np.abs(np.fft.fftshift(np.fft.fft(y-m, n=N, axis=0), axes=0))
    mini = np.argmin(Y, axis=0)
    maxi = np.argmax(Y, axis=0)
    periodAll = N/np.abs(mini-maxi)
    bpm = 60/dt/periodAll
    periodAll[bpm<minbpm] = 0
    periodAll[bpm>maxbpm] = 0
    period = np.amax(periodAll)
    if period == 0:
        print('No waveform: Skipping bad data')
        for j in range(nc):
            data[title + '_waveform'][j]=np.nan*np.zeros((10,))
        return data

    ### separate periods ###
    yp = nc*[0]
    peaks = nc*[0]
    lengths = nc*[0]
    yf = gaussian_filter1d(y, 1, axis=0, mode='reflect')
    for i in range(nc):

        ### get dividing pts
        yp[i] = yf[1:, i]-yf[:-1, i]
        minPeriod = int(np.floor(period-1))
        peaks[i] = find_peaks(yp[i]**2, distance=minPeriod)[0]
        lengths[i]= peaks[i][1:]-peaks[i][:-1]

        ### get segments only pick ones that are close < +/-1 of approximation
        maxlength = np.ceil(period)
        minlength = np.floor(period)
        L = int(maxlength+4)
        usesegment = ( (minlength<=np.array(lengths[i])).astype(int) + (maxlength>=np.array(lengths[i])).astype(int) ) == 2
        if (peaks[i][0]<2):
            usesegment[0]=False
        nsegments =int(np.sum(usesegment))
        if nsegments == 0:
            data[title + '_waveform'][i]=np.nan*np.zeros((10,))
            continue
        segments = np.zeros((L, nsegments))
        ind=0
        order = np.argsort(lengths[i])
        for j in range(len(lengths[i])):
            jj = order[j]
            if lengths[i][jj]==maxlength:
                t1 = peaks[i][jj]-2
                t2 = peaks[i][jj+1]+2
                if t1>=0 and t2<=(y.shape[0]-1):
                    segments[:, ind] = y[t1:t2, i]
                    ind += 1
            elif lengths[i][jj]==minlength:
                t1 = peaks[i][jj]-2
                t2 = peaks[i][jj+1]+3
                if t1>=0 and t2<=(y.shape[0]-1):
                    segments[:, ind] = y[t1:t2, i]
                    ind += 1

        ### subpixel alignment
        moving = segments[2:-2, -1] #align everything to the last one
        ncc = np.zeros((5,))
        alignedsegments = np.zeros((len(moving), nsegments))
        alignedsegments[:, -1] = moving
        xx=np.arange(L)

        for j in range(nsegments-1):
            fixed = segments[:, j]
            for k in range(5):
                out = pearsonr(fixed[k:k+len(moving)], moving)
                ncc[k]=out[0]
            indmax = np.argmax(ncc)
            if indmax>0 and indmax<4:
                a = 0.5*(ncc[indmax-1]-2*ncc[indmax]+ncc[indmax+1])
                b = 0.5*(ncc[indmax+1]-ncc[indmax-1])
                dk = -b/a/2
            else:
                dk = 0
            kmax = indmax + dk
            f = interp1d(xx, fixed, kind='linear')
            xnew = np.arange(kmax, kmax+len(moving))
            if len(xnew)>alignedsegments.shape[0]:
                xnew = xnew[:-1]
            alignedsegments[:, j] = f(xnew)

        data[title + '_waveform'][i] = np.mean(alignedsegments, axis=1)
        data[title + '_waveform_std'][i] = np.std(alignedsegments, axis=1)
        data[title + '_waveform_std_mean'][i] = np.mean(np.std(alignedsegments, axis=1), axis=0)
        data[title + '_waveform_passedCount'][i] = alignedsegments.shape[1] + i/8


    # check that results from the different cameras are aligned
    data[title + '_waveform'], data[title + '_waveform_std'], = aligncameras(data[title + '_waveform'], data[title + '_waveform_std'])

    return data

def auc(data, sign, title):

    nc = len(data[title + '_waveform'])
    x = np.zeros((nc,))
    for j in range(nc):
        f = sign*data[title][j]
        if all(np.isfinite(f)):
            mini = np.amin(f)
            f -= mini
            maxi = np.amax(f)
            f /= maxi
            x[j] = np.mean(f)
    data[title + '_auc'] = x

    return data


def aligncameras(data, data_std):
    # align to the camera that has most in common with the others (nearest pixel)

    m = 5 #how many steps to sample
    nc = len(data)
    ncc = np.zeros((nc, nc))

    # find best camera
    for i in range(nc):
        for j in range(nc):
            bobs = len(data[i])==len(data[j])
            your = all(np.isfinite(data[i]))
            uncle = all(np.isfinite(data[j]))
            if bobs and your and uncle:
                ncc[i, j] = pearsonr(data[i], data[j])[0]
    fixedind = np.argmax(np.nansum(ncc, axis=0))
    fixed = data[fixedind]

    # align the others to it
    for i in range(nc):
        if i != fixedind and all(np.isfinite(data[i])):
            moving = data[i]
            moving_std = data_std[i]
            ncc = np.zeros((2*m,))
            for j in range(-m, m):
                ncc[j+m] = pearsonr(fixed, np.concatenate((moving[j:], moving[:j])))[0]
            ncc[np.isnan(ncc)]=0
            ind = np.argmax(ncc)-m
            data[i] = np.concatenate((moving[ind:], moving[:ind]))
            data_std[i] = np.concatenate((moving_std[ind:], moving_std[:ind]))

    return data, data_std

def amplitude(d, title):

    nc = len(d[title + '_mean'])
    d['amplitude']=np.zeros((nc,))
    d[title + '_waveform_AmpOverStd_mean'] = np.zeros((nc,))
    for i in range(nc):
        d['amplitude'][i] = np.amax(d[title + '_waveform'][i]) - np.amin(d[title + '_waveform'][i])
    d['modulationdepth'] = d['amplitude']/d[title + '_mean']
    d[title + '_waveform_AmpOverStd_mean'] = d['amplitude'] / d[title + '_waveform_std_mean']

    return d

def ratio(data, title):

    data[title + '_ratio'] = data[title][-1]/data[title][0]

    return data

def combineAndPlotPositions(data, data_dict, title, batch):

    npositions = len(data)
    nc = len(data[0][title])
    normdiff, diff, d = getdifference(data, title, data_dict['col'],  data_dict['row'])
    normpiff, piff, p = getdifference(data, 'pulse', data_dict['col'],  data_dict['row'])
    diff *= np.isfinite(piff)

    x = combinepositions(data, title)
    y = combinepositions(data, 'pulse')
    x *= np.isfinite(y)

    if data_dict['scanorder'] == 'hhc30':

        forehead_left = np.nanmean(x[0:6, :])
        fissure_left  = np.nanmean(x[6:9, :])
        temple_left   = np.nanmean(x[9:15, :])
        forehead_right = np.nanmean(x[15:21, :])
        fissure_right  = np.nanmean(x[21:24, :])
        temple_right   = np.nanmean(x[24:30, :])

        data[0][title + '_combined'] = np.array([forehead_left, fissure_left, temple_left, forehead_right, fissure_right, temple_right])

        forehead = diff[0:6, :]
        forehead_camera = np.nanmean(forehead, axis=0)
        forehead_position = np.nanmean(forehead, axis=1)
        forehead_all = np.nanmean(forehead)

        fissure = diff[6:9, :]
        fissure_camera = np.nanmean(fissure, axis=0)
        fissure_position = np.nanmean(fissure, axis=1)
        fissure_all = np.nanmean(fissure)

        temple = diff[9:, :]
        temple_camera = np.nanmean(temple, axis=0)
        temple_position = np.nanmean(temple, axis=1)
        temple_all = np.nanmean(temple)

        data[0][title + '_combined_diff'] = np.array([forehead_all, fissure_all, temple_all])

        if batch == False:
            fig, ax = plt.subplots(nrows=2, ncols=3)
            fig.suptitle(title)
            ax[0, 0].plot(data_dict.get('rho'), forehead_camera, '.-')
            ax[0, 0].plot(data_dict.get('rho'), forehead_all*np.ones((nc,)), '--')
            ax[0, 0].set_ylabel('Left - Right')
            ax[0, 0].set_xlabel('SD separation (mm)')
            ax[0, 0].set_title('Forehead')
            ax[1, 0].plot(forehead_position, '.-')
            ax[1, 0].plot(forehead_all*np.ones((len(forehead_position),)), '--')
            ax[1, 0].set_ylabel('Left - Right')
            ax[1, 0].set_xlabel('Positions')
            ax[0, 1].plot(data_dict.get('rho'), fissure_camera, '.-')
            ax[0, 1].plot(data_dict.get('rho'), fissure_all*np.ones((nc,)), '--')
            ax[0, 1].set_xlabel('SD separation (mm)')
            ax[0, 1].set_title('Fissure')
            ax[1, 1].plot(fissure_position, '.-')
            ax[1, 1].plot(fissure_all*np.ones((len(fissure_position),)), '--')
            ax[1, 1].set_xlabel('Positions')
            ax[0, 2].plot(data_dict.get('rho'), temple_camera, '.-')
            ax[0, 2].plot(data_dict.get('rho'), temple_all*np.ones((nc,)), '--')
            ax[0, 2].set_xlabel('SD separation (mm)')
            ax[0, 2].set_title('Temple')
            ax[1, 2].plot(temple_position, '.-')
            ax[1, 2].plot(temple_all*np.ones((len(temple_position),)), '--')
            ax[1, 2].set_xlabel('Positions')

        return data

    else:
        print('Nothing programmed to combine this scan order')

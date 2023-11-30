# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 08:00:15 2020

@author: soren
"""

import csv
import numpy as np
#import re
import matplotlib.pyplot as plt
import copy
from scipy.ndimage import gaussian_filter1d
from scipy.signal import medfilt
from scipy.optimize import least_squares


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

def removedark2(data):
    
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
    
    return d

def croptime_proto(d, data_dict):

    nstart = data_dict.get('nstart', 0)
    nstop = data_dict.get('nstop', d['mean'].shape[0]+1)
    for key in d.keys():
        d[key] = d[key][nstart:nstop, :]

    return d

def shotnoise(d, data_dict):
    
    d['sigma'] = np.sqrt( d['sigma']**2 - data_dict['K']*data_dict['gain']*d['mean'] )
    d['contrast'] = d['sigma']/d['mean']
    
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
                if d[key][i]==0:
                    d[key][i]=np.nan
                
    return d            
    
def plotscandataintime(d):
    
    mean = np.zeros((1200, 4))
    std = np.zeros((1200, 4))
    contrast = np.zeros((1200, 4))
    for i in range(12):
        m = int(i*100)
        n = int((i+1)*100)
        mean[m:n, :] = d[i]['mean']
        std[m:n, :] = d[i]['std']
        contrast[m:n, :] = d[i]['contrast']
    
    c = 'rgbcmkkmcbgr'
    m = '<<<<<<>>>>>>' 

    fig, ax = plt.subplots(nrows=4, ncols=3)
    for j in range(12):
        for i in range(4):
            ax[i, 0].plot(np.arange(j*100, (j+1)*100), d[j]['mean'][:, i], c[j]+m[j])
            ax[i, 1].plot(np.arange(j*100, (j+1)*100), d[j]['std'][:, i], c[j]+m[j])
            ax[i, 2].plot(np.arange(j*100, (j+1)*100), d[j]['contrast'][:, i], c[j]+m[j])
            ax[i, 0].plot(j*np.array([100, 100]), np.array([np.amin(mean[:, i], axis=0), np.amax(mean[:, i], axis=0)]), ':r')
            ax[i, 1].plot(j*np.array([100, 100]), np.array([np.amin(std[:, i], axis=0), np.amax(std[:, i], axis=0)]), ':r')
            ax[i, 2].plot(j*np.array([100, 100]), np.array([np.amin(contrast[:, i], axis=0), np.amax(contrast[:, i], axis=0)]), ':r')
    fig.suptitle('Mean              Std              Contrast')

    fig, ax = plt.subplots(nrows=4, ncols=2)
    for j in range(12):
        for i in range(4):
            ax[i, 0].plot(d[j]['mean'][:, i], d[j]['std'][:, i],      c[j]+m[j])
            ax[i, 1].plot(d[j]['mean'][:, i], d[j]['contrast'][:, i], c[j]+m[j])
            ax[i, 0].set_xlabel('Mean')
            ax[i, 0].set_ylabel('Std')
            ax[i, 1].set_xlabel('Mean')
            ax[i, 1].set_ylabel('Contrast')
            
    return mean, std, contrast        
         

def getoffsets(d):
    
    fig, ax = plt.subplots(nrows=1, ncols=4)
    for i in range(2):
        for j in range(4):
            ax[j].plot((d['std'][:, 0]-i)/d['mean'][:, j])
    


def getpulse(d, data_dict, dt):
    
    N = 512
    m = np.mean(d['contrast'], axis = 0)
    d['pcontrast'] = np.abs(np.fft.fftshift(np.fft.fft(d['contrast']-m, n=N, axis=0), axes=0))
    #N = d['contrast'].shape[0]
    p = d['pcontrast'][int(N/2):, np.argmin(data_dict['rho'])]  #closest camera 
    ind = np.argmax(p)
    if ind>20:
        d['pulse'] = 60*ind/N/dt
    else:
        d['pulse'] = np.nan
    
    
    return d
       
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

def fitattenuation(mean, p, pp, musp, sigma):
    
    nt = mean.shape[0]
    nexp = mean.shape[2]

    logdata = np.log((p**2)*mean[0, :, 0])
    M = M = np.stack((np.ones(len(p),), -p), axis=1)
    Minv = np.linalg.pinv(M)
    temp = Minv @ logdata
    initial_guess = np.array([np.exp(temp[0]), temp[1]])
    #print(initial_guess)

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
                mean_calc[i, :, j] = I0[i, j]*G_semiinf(mu_eff[i, j], pp) #*norm
            else:
                mu_eff[i, j] = np.nan
                I0[i, j] = np.nan
                mua[i, j] = np.nan
                mean_calc[i, :, j] = np.nan

    colors = 'brgcmykbrgcmykbrgcmykbrgcmyk'
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
    plt.title('u_eff = ' + np.array2string(np.squeeze(mu_eff), precision = 2) + ' mm-1')#,  u_a = ' +  np.array2string(np.mean(mua), precision = 2) + ' mm-1')
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
    plt.title('u_eff = ' + np.array2string(np.squeeze(mu_eff), precision = 2) + ' mm-1')#,  u_a = ' +  np.array2string(np.mean(mua), precision = 2) + ' mm-1')
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
        
    colors = 'brgcmykbrgcmykbrgcmykbrgcmyk'
    num1 = 0
    num2 = 0
    nexp = len(params_dict['exp'])
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
    plt.title('Blood Flow Index = ' + np.array2string(alphaDb, precision = 2) + ' mm2/s')
    plt.legend(('fit', 'data'))

    return alphaDb, beta, offset

def residual_speckle_multi(x, rhoplus, data):
    
    exps = rhoplus[-data.shape[1]:]
    expinds = (np.round(exps/rhoplus[5])).astype('int32')
    rhoplus = rhoplus[:-data.shape[1]]
    calculated = specklecontrast_numerical_fitAlphaDb(x, rhoplus)
    R = calculated[expinds, :] - data.T
    
    return R.flatten()
    

def plotscandata(data, data_dict, title):

    row = data_dict.get('row')
    col = data_dict.get('col')    
    ns = len(data)
    fig, ax = plt.subplots(nrows=9, ncols=int(ns/2))
    for i in range(ns):
        ax[int(row[i]*5+0), int(col[i])].plot(data[i][title][:, 0], '.-')
        ax[int(row[i]*5+1), int(col[i])].plot(data[i][title][:, 1], '.-')
        ax[int(row[i]*5+2), int(col[i])].plot(data[i][title][:, 2], '.-')
        ax[int(row[i]*5+3), int(col[i])].plot(data[i][title][:, 3], '.-')
    fig.suptitle(title)
    for i in range(int(ns/2)):
        ax[4, i].axis('off')

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
    
    normdiff = np.zeros((hh, 4))
    left = np.arange(hh)
    right = hh+np.flip(left)
    for i in range(hh):
        normdiff[i, :] = 2*(data[left[i]][title]-data[right[i]][title])/(data[left[i]][title]+data[right[i]][title])
    
    fig, ax = plt.subplots(nrows=2, ncols=1)
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
    ax[0].set_title(title)
    ax[0].set_ylabel('Back of Head')
    ax[0].set_xlabel('Right Side of Head')

    ax[1].plot(normdiff, '.-')
    ax[1].legend(data_dict['rho'])
    ax[1].set_title('Left Right Difference')
    ax[1].set_xlabel('<-- Back            Front -->')
    ax[1].set_ylabel('2 * (L-R)/(L+R)')       
    ax[1].set_xticks([])
    
    
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
    
    normdiff = np.zeros((hh,))
    left = np.arange(hh)
    right = hh+np.flip(left)
    for i in range(hh):
        normdiff[i] = 2*(data[left[i]][title]-data[right[i]][title])/(data[left[i]][title]+data[right[i]][title])
    
    fig, ax = plt.subplots(nrows=2, ncols=1)
    im = np.nan*np.zeros((Ny, Nx))
    for i in range(ns):
        xx = int(x[i])
        yy = int(y[i])
        im[yy-h:yy+h, xx-h:xx+h] = data[i][title]
        if np.isfinite(data[i][title]):
            value = np.array2string(data[i][title], precision=2)
        else:
            value = 'nan'
        ax[0].text(xx-hh*1.5, yy, value, color=[1, 1, 1])
        
    ax[0].imshow(im)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_title(title)
    ax[0].set_ylabel('Back of Head')
    ax[0].set_xlabel('Right Side of Head')

    ax[1].plot(normdiff, '.-')
    #ax[1].legend(data_dict['rho'])
    ax[1].set_title('Left Right Difference')
    ax[1].set_xlabel('<-- Back            Front -->')
    ax[1].set_ylabel('2 * (L-R)/(L+R)')       
    ax[1].set_xticks([])    
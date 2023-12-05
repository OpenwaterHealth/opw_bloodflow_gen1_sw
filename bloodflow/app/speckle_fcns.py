# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 20:31:05 2020

@author: soren
"""

import numpy as np
from scipy.optimize import least_squares

def get_rhoplus_multi(params_dict):
  ''' create rho+ data vector from params_dict '''
  nrho = len(params_dict['rho'])
  nexp = len(params_dict['exp'])

  rhoplus = np.zeros((8+nrho+nexp,))
  rhoplus[0] = params_dict['n']    # index of refraction
  rhoplus[1] = params_dict['wv']   # wavelength, mm
  rhoplus[2] = params_dict['mua']  # absorption coefficient, mm^-1
  rhoplus[3] = params_dict['musp'] # reduced scattering coefficient, mm^-1
  rhoplus[4] = params_dict['Tmax'] # maximum exposure time being simulated, s
  rhoplus[5] = params_dict['dtau'] # discretization of time step length, s
  rhoplus[6] = params_dict['zb']   # extrapolated boundary (for method of images), mm
  rhoplus[7] = params_dict['beta'] # constant from 0 to 1 depending on pixel size, polarization,
                                   # and abberations, unitless
  rhoplus[8:(8+nrho)] = params_dict['rho'] # distance between source and detector fibers, mm
  rhoplus[(8+nrho):] = params_dict['exp']

  return rhoplus

def specklecontrast_numerical_fitAlphaDb(x, rhoplus):
  ''' fit speckle contrast numerically from:
  args:
      x:
      rhoplus: augmented data vector
  '''
  n = rhoplus[0]    # index of refraction
  wv = rhoplus[1]   # wavelength, mm
  mua = rhoplus[2]  # absorption coefficient, mm^-1
  musp = rhoplus[3] # reduced scattering coefficient, mm^-1
  Tmax = rhoplus[4] # maximum exposure time being simulated, s
  dtau = rhoplus[5] # discretization of time step length, s
  zb = rhoplus[6]   # extrapolated boundary (for method of images), mm
  alpha = 1         # fraction of scatterers that are moving
  if len(x) == 1:
    beta = rhoplus[7]
    Db = x[0]     # diffusion coefficient of the scatterers, mm^2/s
    offset = 0
  elif len(x) == 2:
    beta = x[0]   # constant from 0 to 1 depending on pixel size, polarization, and abberations,
                  # unitless
    Db = x[1]
    offset = 0
  elif len(x) == 3:
    beta = x[0]
    Db = x[1]
    offset = x[2]
  rho = rhoplus[8:] # distance between source and detector fibers, mm

  ### calculations ###
  k0 = 2*np.pi*n/wv      #wavenumber
  tau = np.arange(0, Tmax+dtau, dtau)
  dr2 = (0.5*(tau[1:]+tau[:-1]))*(6*Db)
  z0 = 1/musp
  contrast = np.zeros((len(dr2), len(rho)))   # (exposures, rhos)
  r1 = np.sqrt(rho**2 + z0**2)                # len = rhos
  r2 = np.sqrt(rho**2 + (z0+(2*zb))**2)       # len = rhos
  k = np.sqrt(alpha*(k0**2)*(musp**2)*dr2 + 3*mua*musp) # len = exposures
  sourceterm = np.exp(-np.outer(k, r1))/np.tile(r1, (len(dr2), 1))
  imageterm = np.exp(-np.outer(k, r2))/np.tile(r2, (len(dr2), 1))
  G1 = sourceterm - imageterm
  g1 = np.zeros(G1.shape)  # (exposures, rhos)
  for i in range(len(rho)):
    g1[:, i] = G1[:, i] / G1[0, i]
  taus = np.tile(tau[:, np.newaxis], (1, len(rho)))
  part1 = np.cumsum((g1**2), axis=0) / taus[1:]
  part2 = np.cumsum(0.5*(taus[1:]+taus[:-1])*(g1**2), axis=0) / (taus[1:]**2)
  part1minuspart2 = part1 - part2
  contrast = np.sqrt(2 * beta * dtau * part1minuspart2) + offset

  return contrast

def residual_speckle_multi(x, rhoplus, data):
  ''' calculate residual speckle from:
  args:
      x: minimization arg from least_squares
      rhoplus: augmented parameter vector
  '''
  exps = rhoplus[-data.shape[1]:]
  expinds = (np.round(exps/rhoplus[5])).astype('int32')
  rhoplus = rhoplus[:-data.shape[1]]
  calculated = specklecontrast_numerical_fitAlphaDb(x, rhoplus)
  R = calculated[expinds, :] - data.T

  return R.flatten()

def G_semiinf(k, rho):
  ''' return G semi-infinite function of:
  args:
      k: single parameter value from least_squares
      rho: source-to-sensor distances
  '''
  G = np.exp(-k*rho)/(rho**2)

  return G

def residual_intensity2(x, rho, data):
  ''' calculate residual intensity from:
  args:
      x: minimization arg from least_squares
      rho: distances between source and detectors
      data: data range to fit
  '''
  if len(x) == 2:
    R = x[0]*G_semiinf(x[1], rho) - data
  elif len(x) == 1:
    R = G_semiinf(x[0], rho) - data
  return R

def fitattenuation(mean, p, nt, nexp, musp, initial_guess):
  ''' find mu_eff,mua attenuation from:
  args:
      mean: data mean intensity
      p: distances between source and detectors
      nt: number of time points
      nexp: number of exposures
      musp: reduced scattering coefficient, mm^-1
      initial_guess: starting point for solution
  '''
  I0 = np.zeros((nt, nexp))
  fitI0 = (len(initial_guess) == 2)

  mu_eff = np.zeros((nt, nexp))
  mua = np.zeros((nt, nexp))

  for i in range(nt):
    for j in range(nexp):
      data = np.squeeze(mean[i, :, j])
      out = least_squares(
        residual_intensity2, initial_guess, bounds=(0, np.inf), args=(p, data))
      if fitI0:
        mu_eff[i, j] = out.x[1]
        I0[i, j] = out.x[0]
      else:
        mu_eff[i, j] = out.x[0]
        I0[i, j] = 1
      mua[i, j] = mu_eff[i, j]**2 / (3*musp)

  return mu_eff, mua

def fitflow_multi(params_dict, contrast, nt, mua):
  ''' fit the flow rate from:
  args:
      param_dict: static parameters
      contrast: speckle contrast values
      nt: number of time points
      mua: as returned from fitattenuation
  '''
  if params_dict['fitbeta'] == 2:
    initial_guess = [0.3, 1e-6, 0]    #[beta, alpha*Db]
    lb = [0, 0, 0]
    ub = [1, np.inf, np.inf]
  elif params_dict['fitbeta'] == 1:
    initial_guess = [0.3, 1e-6]    #[beta, alpha*Db]
    lb = [0, 0]
    ub = [1, np.inf]
  else:
    initial_guess = [1e-6]
    lb = 0
    ub = np.inf

  alphaDb = np.zeros((nt, 1))
  beta = np.zeros((nt, 1))
  offset = np.zeros((nt, 1))
  rhoplus = get_rhoplus_multi(params_dict)

  for i in range(nt):
    for j in range(alphaDb.shape[1]):  # pylint: disable=unsubscriptable-object
      if params_dict['usemueff'] == 1:
        params_dict['mua'] = np.mean(mua[i, :])
        rhoplus = get_rhoplus_multi(params_dict)
      data = contrast[i, :, :]
      out = least_squares(
        residual_speckle_multi, initial_guess, bounds=(lb, ub), args=(rhoplus, data))

      if params_dict['fitbeta'] == 2:
        beta[i, j] = out.x[0]
        alphaDb[i, j] = out.x[1]
        offset[i, j] = out.x[2]
      if params_dict['fitbeta'] == 1:
        beta[i, j] = out.x[0]
        alphaDb[i, j] = out.x[1]
      else:
        beta[i, j] = params_dict['beta']
        alphaDb[i, j] = out.x[0]

  return alphaDb, beta, offset

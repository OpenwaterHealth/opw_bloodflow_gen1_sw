''' unit test for speckle_fcns '''

import numpy as np
import speckle_fcns as fcns

# example data
mean = np.array([0.083, 0.4827, 2.9799, 45.5313])     # mean intensity
contrast = np.array([0.1206, 0.1406, 0.1933, 0.3021]) # speckle contrast
rhos = np.array([38, 28.5, 19, 9.5]) # distances between source and detectors
exps = np.array([200])*1e-6          # pulse widths / exposures
nt = 1                               # number of time points

# parameters
params = {
  'n': 1.4,      # index of refraction
  'wv': 8.5e-4,  # wavelength, mm
  'mua': 0.005,  # absorption coefficient, mm^-1
  'musp': 1,     # reduced scattering coefficient, mm^-1
  'Db': 1e-6,    # alpha * diffusion coefficient of the scatterers, mm^2/s
  'dV2': 1,      # mean squared velocity of the scatterers, (mm/s)^2
  'Tmax': 1e-3,  # maximum exposure time being simulated, s
  'dtau': 1e-7,  # discretization of time step length, s
  'rho': rhos,   # distance between source and detector fibers, mm
  'exp': exps,   # exposure times, s
  'zb': 0.1,     # extrapolated boundary (for method of images), mm
  'beta': 0.25,  # constant from 0 to 1 depending on pixel size, polarization, and abberations, unitless
  'initial_guess': [40, 0.1], # for fitting I0 and mu_eff
  'usemueff': True,     # use the fitted mueff value in the flow fitting
  'fitbeta': 0,         # 1=fit for the beta parameter, 2=fit for beta and an offset
}

# arrange data in form expected by fitting functions (time, distance, pulse width)
mean = np.expand_dims(np.expand_dims(mean, axis=1), axis=0)
contrast = np.expand_dims(np.expand_dims(contrast, axis=1), axis=0)

p = np.squeeze(params['rho'])
exp = params['exp']
eps = 1e-10

def test_fitattenuation():
  mu_eff, mua = fcns.fitattenuation(mean, p, nt, len(exp), params['musp'], params['initial_guess'])
  assert abs(mu_eff[0][0] - 0.13995480376532396) < eps
  assert abs(mua[0][0] - 0.006529115698996779) < eps

def test_fitflow_multi():
  mu_eff, mua = fcns.fitattenuation(mean, p, nt, len(exp), params['musp'], params['initial_guess'])
  alphaDb, beta, offset = fcns.fitflow_multi(params, contrast, nt, mua)
  assert abs(alphaDb[0][0] - 1.1254315142302882e-06) < eps
  assert abs(beta[0][0] - 0.25) < eps
  assert offset == 0

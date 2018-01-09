#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  seismogram_noise.py
#   Purpose:   Entry point into seismogram_noise
#   Author:    Simon Staehler
#   Email:     mail@simonstaehler.com
#   License:   MIT License
# -------------------------------------------------------------------

import numpy as np
from scipy.fftpack import fft, ifft, fftfreq
import sys

def create_noise(dt, npts, f_in, power_in, 
                 interpolate='loglog'):
    """
    create noise series with given power spectrum
    
    Keywords:
    :type  dt: float
    :param dt: desired dt in seconds

    :type  npts: int
    :param npts: desired number of samples
    
    :type  f_in: numpy.array
    :param f_in: frequency array of input power spectrum

    :type  power_in: numpy.array
    :param power_in: input power spectrum

    :type  interpolate: 'linear' or 'loglog'
    :param interpolate: interpolation can either be done linear or
                        in a log-log plot. The latter is usually 
                        preferred, since noise is assumed to be 
                        linear between points in log-log plot.
    """

    if ((power_in<0).any()):
        raise ValueError('Power spectrum must be >=0 everywhere')

    if ((f_in<0).any()):
        raise ValueError('Frequency vector must be >=0 everywhere')
    
    # Create random time series of desired length
    noise = np.random.randn(npts)
    noise_fd = fft(noise)
    
    # Normalize noise spectrum
    noise_fd /= np.abs(noise_fd)
    
    # Get frequency vector
    f = fftfreq(npts, d=dt) 

    # calculate energy spectral density
    energy_in = np.sqrt(power_in)
    
    # Multiply with desired spectrum
    if interpolate == 'linear':
        noise_amp_ipl = np.interp(x=abs(f), 
                                  xp=f_in, 
                                  fp=energy_in)
    elif interpolate == 'loglog':
        # Remove exact zeros from energy spectrum
        energy_in[energy_in==0] = sys.float_info.min

        # Since the self-noise is usually defined as linear in 
        # log-log plots, we need to interpolate the logarithms
        noise_amp_ipl = 10**np.interp(x=np.log10(abs(f)), 
                                      xp=np.log10(f_in), 
                                      fp=np.log10(energy_in))
    else: 
        raise ValueError('Unknown interpolation scheme %s' % interpolate)
    noise_fd *= np.abs(noise_amp_ipl)
    noise_td = ifft(noise_fd * np.sqrt(npts / 2 / dt))
    
    return noise_td.real

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
import os
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm


def external_models():
    data_path = os.path.join(os.path.dirname(__file__),
                             'data')
    fnams = dict()
    fnams['Tcompact'] = os.path.join(data_path,
                                     'noise_Trillium_compact.txt')
    fnams['STS2'] = os.path.join(data_path,
                                 'noise_STS2.txt')
    fnams['CMG40T-OBS'] = os.path.join(data_path,
                                       'noise_CMG40T_OBS.txt')
    fnams['external'] = None
    fnams['NHNM'] = 'routine'
    fnams['NLNM'] = 'routine'

    return fnams


def get_spectrum(model):
    fnams = external_models()
    if fnams[model]:
        if model == 'NHNM':
            p, power = get_nhnm()
            f_in = 1./p
            power_in = 10**(power/10)

        elif model == 'NLNM':
            p, power = get_nlnm()
            f_in = 1./p
            power_in = 10**(power/10)

        else:
            spec = np.loadtxt(fnams[model])
            f_in = spec[:, 0]
            power_in = spec[:, 1]**2

    return f_in, power_in


def add_noise(st, model='external', f_in=None, power_in=None, **kwargs):
    """
    add noise series with given power spectrum to Stream object

    Keywords:
    :type  st: obspy.Stream
    :param dt: ObsPy Stream object where noise should be added

    :type  model: string
    :param model: Noise model name. Allowed options are
                  Tcompact:   Trillium compact
                  STS2:       Streckeisen STS-2
                  CMG40T-OBS: OBS version of Guralp CMG-40T
                  NHNM:       Petersen New High Noise model
                  NLNM:       Petersen New Low Noise model
                  external:   Own noise model provided by variables
                              f_in and power_in

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
    if model == 'external':
        if f_in is None or power_in is None:
            raise ValueError('Either specify a noise model or provide one')
    else:
        f_in, power_in = get_spectrum(model)

    for tr in st:
        tr.data += create_noise(dt=tr.stats.delta,
                                npts=tr.stats.npts,
                                f_in=f_in,
                                power_in=power_in,
                                **kwargs)

    return st


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

    if ((power_in < 0).any()):
        raise ValueError('Power spectrum must be >=0 everywhere')

    if ((f_in < 0).any()):
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
    noise_amp_ipl = interpolation(f, f_in, energy_in, interpolate)
    noise_amp_ipl[f == 0] = np.interp(x=0.0,
                                      xp=f_in,
                                      fp=energy_in)

    noise_fd *= np.abs(noise_amp_ipl)
    noise_td = ifft(noise_fd * np.sqrt(npts / 2 / dt))

    return noise_td.real


def interpolation(f, f_in, energy_in, interpolate='loglog'):
    if interpolate == 'linear':
        noise_amp_ipl = np.interp(x=abs(f),
                                  xp=f_in,
                                  fp=energy_in)
    elif interpolate == 'loglog':
        # Remove exact zeros from energy spectrum
        energy_in[energy_in == 0] = sys.float_info.min

        # Since the self-noise is usually defined as linear in
        # log-log plots, we need to interpolate the logarithms
        noise_amp_ipl = 10**np.interp(x=np.log10(abs(f)),
                                      xp=np.log10(f_in),
                                      fp=np.log10(energy_in))
    else:
        raise ValueError('Unknown interpolation scheme %s' % interpolate)

    return noise_amp_ipl


def get_noise_power(freq, model='external', interpolate='loglog',
                    f_in=None, power_in=None):
    if model == 'external':
        if f_in is None or power_in is None:
            raise ValueError('Either specify a noise model or provide one')
    else:
        f_in, power_in = get_spectrum(model)

    energy_in = np.sqrt(power_in)

    return interpolation(freq, f_in, energy_in, interpolate)

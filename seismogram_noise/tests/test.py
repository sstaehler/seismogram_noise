import obspy
import seismogram_noise
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np


for model in ['STS2', 'Tcompact', 'CMG40T-OBS', 'NHNM', 'NLNM']:
    # Read data (basically zeros)
    st = obspy.read('../data/testdata.mseed')

    # Add noise
    seismogram_noise.add_noise(st, model=model)
    power, f = mlab.psd(x=st[0].data, Fs=st[0].stats.sampling_rate, NFFT=2**14)
    plt.plot(1./f[f>0], 10*np.log10(power[f>0]), label='resulting spectrum')

    # Add theoretical spectrum
    spec = seismogram_noise.get_spectrum(model)
    plt.plot(1./spec[0], 10*np.log10(spec[1]), 'r', label='theory', lw=2)

    plt.xlim(0.1, 500)
    plt.xscale('log')
    plt.xlabel('period / seconds')
    plt.ylabel('amplitude / dB')
    plt.legend()
    plt.savefig('testspec_%s.png' % model)
    plt.close()

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.gridspec as gridspec
from PyAstronomy.pyasl import binningx0dt

def combi_plotter(lc, filter=None, save_pth=None):
    
    ti= lc[1].data['time (BJD-TDB)'] #-2450000
    fl= lc[1].data['flux']
    fl_err= lc[1].data['sflux']
    
    fwhm= lc[1].data['fwhm']
    peak_flux= lc[1].data['peak']
    air_mass= lc[1].data['airmass']
    sky= lc[1].data['bkg']
    xshift= lc[1].data['dx']
    yshift= lc[1].data['dy']
    exptime= lc[1].data['exptime']
    defoc= lc[1].data['defoc']

    fig = plt.figure(figsize = (8, 12))
    gs = gridspec.GridSpec(4,2, width_ratios = [1, 1], height_ratios = [1, 1, 1, 1])
    
    ax0=plt.subplot(gs[0,0])
    ax0.errorbar(ti, fl, yerr=fl_err, fmt='o', c='k', ms = 2)
    ax0.set_xlabel('Time [BJD-TDB]', fontsize=10)
    ax0.set_ylabel('Relative flux', fontsize=10)
    
    
    ax1=plt.subplot(gs[0, 1])
    ax1.scatter(ti, peak_flux, s=2.5, c= 'navy')
    ax1.set_xlabel('Time [BJD-TDB]', fontsize=10)
    ax1.set_ylabel('Peak flux', fontsize=10)
    
    
    ax2=plt.subplot(gs[1,0])
    ax2.scatter(ti, air_mass, s=2.5, c= 'navy')
    ax2.set_xlabel('Time [BJD-TDB]', fontsize=10)
    ax2.set_ylabel('Air mass', fontsize=10)
    
    ax3=plt.subplot(gs[1,1])
    ax3.scatter(ti, sky, s=2.5, c= 'navy')
    ax3.set_xlabel('Time [BJD-TDB]', fontsize=10)
    ax3.set_ylabel('Sky', fontsize=10)
    
    ax4= plt.subplot(gs[2,0])
    ax4.scatter(ti, xshift, s=2.5, c= 'green')
    ax4.scatter(ti, yshift, s=2.5, c= 'red')
    ax4.set_xlabel('Time [BJD-TDB]', fontsize=10)
    ax4.set_ylabel('Shifts [pixel]', fontsize=10)
    
    ax5= plt.subplot(gs[2,1])
    ax5.scatter(ti, fwhm, s= 2.5, c= 'navy')
    ax5.set_xlabel('Time [BJD-TDB]', fontsize=10)
    ax5.set_ylabel('FWHM of Gaussian fit [pixels]', fontsize=10)
    
    ax6= plt.subplot(gs[3,0])
    ax6.scatter(ti, exptime, s=2.5, c= 'navy')
    ax6.set_xlabel('Time [BJD-TDB]', fontsize=10)
    ax6.set_ylabel('Exposure time', fontsize=10)
    
    ax7= plt.subplot(gs[3,1])
    ax7.scatter(ti, defoc, s=2.5, c= 'navy')
    ax7.set_xlabel('Time [BJD-TDB]', fontsize=10)
    ax7.set_ylabel('Defocus [mm]', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(save_pth+f'combi_{filter}.png', dpi=200)
    plt.show()
    
    if save_pth:
        np.savetxt(save_pth+f'lc_data_{filter}.dat', np.transpose([ti, fl, fl_err, xshift, yshift, air_mass, fwhm, sky, exptime, peak_flux]), fmt= '%3.5f', header='time  flux    fluxerr xshift  yshift  airmass fwhm    sky exptime maxval')
    
        
    return
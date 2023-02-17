# Functions used in explore.ipynb
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from scipy.ndimage.filters import gaussian_filter
from astropy.cosmology import FlatLambdaCDM
from copy import copy
from numpy import nan

from pylab import *
import astropy.io.fits as fits
from scipy.io import readsav
from matplotlib.colors import LogNorm
from astropy.visualization import (ZScaleInterval, ImageNormalize)
import matplotlib.patches as patches
#from plot_grids import *
import scipy.integrate as integrate

Red = '\033[91m'
Green = '\033[92m'
Blue = '\033[94m'
Cyan = '\033[96m'
Yellow = '\033[93m'
Magenta = '\033[95m'
CEND = '\033[0m'

# ions wavelenghts
o2 = [3726.03, 3728.82]
mg2 = [2795.5301, 2802.7056]
o3 = [4958.911,5006.843]

#Redshift of the source
z = 0.043118

#check that the python library gives the same result as the  manual calculation
cosmo = FlatLambdaCDM(H0=67.4, Om0=0.315, Tcmb0=2.725)
arcsec_kpc = cosmo.arcsec_per_kpc_proper(z) 
deg_kpc = arcsec_kpc / 60. /60. # degree size of 1 kpc

print(Green + 'Check using FlatLambdaCDM, scale:' + CEND, arcsec_kpc, 1/arcsec_kpc)

#Set some constants
h0 = 67.4
Om = 0.315
Ol = 1 - Om
c_light = 299792.458

#Estimate the Angular Diameter Distance 
def func(x):
    return 1/np.sqrt( Om * (1+x)**3 + Ol)

result = integrate.quad(lambda x: func(x), 0, z)
da = c_light/h0*1/(1+z)* (result[0])
print(Red + 'Angular Diameter Distance is' + CEND, da, Red + 'Mpc' + CEND)

# 4.84814e-6 is rad in arcsec, 1000 is to have the scale in kpc
x = da * 1000 * 4.84814e-6

# 8.096e-5 is the pixel scale along slice, 3600 is to convert it in arcsec (0.291456)
pixel_size = 8.096e-5 * 3600
pixel_scale = x * pixel_size


arcsec_pixel =  pixel_scale / x

####################################################
# Functions to get data for maps

def set_sn(arr, arr2):
    sn_cut = np.where(arr < 3) #change here the SN threshold, now set to 3
    arr2[sn_cut] = nan   #if the SN is less than 5, remove that spaxel
    return arr2

def set_bad(arr):
    bad = (arr > 1e90)   
    arr[bad] = nan   
    return arr   

def get_flux_uncorr(idx, arr, arr2):
    f = copy(arr[idx][1])
    f_err = copy(arr2[idx][1])
    f = set_bad(f)
    f_err = set_bad(f_err)  
    for row in range(f_err.shape[0]):
        for column in range(f_err.shape[1]):
            if(f_err[row][column] < err_thresh):
                f_err[row][column] = err_thresh
    sn = f/f_err
    f = set_sn(sn,f)  
    return f

def get_err_flux(idx, arr):
    f_err = copy(arr[idx][1])
    f_err = set_bad(f_err)  
    return f_err

    
def get_variable(idx, arr, arr2, arr3):
    f = copy(arr[idx][1])
    f_err = copy(arr2[idx][1])
    f = set_bad(f)
    f_err = set_bad(f_err)
    for row in range(f_err.shape[0]):
        for column in range(f_err.shape[1]):
            if(f_err[row][column] < err_thresh):
                f_err[row][column] = err_thresh
    
    sn = f/f_err  
    v = copy(arr3[idx][1])
    v = set_bad(v)
    v = set_sn(sn,v)
    return v

def get_sn(idx, arr, arr2):
    f = copy(arr[idx][1])
    f_err = copy(arr2[idx][1])
    f = set_bad(f)
    f_err = set_bad(f_err)
    for row in range(f_err.shape[0]):
        for column in range(f_err.shape[1]):
            if(f_err[row][column] < err_thresh):
                f_err[row][column] = err_thresh
    
    sn = f/f_err  
    sn_cut = np.where(sn < 3) #change here the SN threshold, now set to 3
    sn[sn_cut] = nan 
    return sn



####################################################
# Calculate error threshold

def err_map(idx, arr, title = None, kpc=False, arcsec=False):
    figure(figsize=(4,4))
    ax= subplot(111)
    o2_err = get_err_flux(idx, arr) 

    im = ax.imshow(o2_err, origin='lower', cmap='viridis', norm=LogNorm())
    colorbar(im)
    ax.title.set_text(title)
    
    x_tick = (17, 37, 57)
    ax.set_xticks(x_tick)
    y_tick = (27, 47, 67)
    ax.set_yticks(y_tick)    
    
    if kpc == True:
        x_labels = np.round(np.array([-20, 0 , 20])*pixel_scale,1)
        ax.set_xticklabels(x_labels,size=10)
        y_labels = np.round(np.array([-20, 0 , 20])*pixel_scale,1)
        ax.set_yticklabels(y_labels,size=10)
        #plt.setp(ax, xlim=custom_xlim, ylim=custom_ylim)
        ax.set_xlabel(r'kpc',fontsize=11)
        ax.set_ylabel(r'kpc',fontsize=11)
    if arcsec == True:
        x_labels = np.round(np.array([-20, 0 , 20])*arcsec_pixel,1)
        ax.set_xticklabels(x_labels,size=10)
        y_labels= np.round(np.array([-20, 0 , 20])*arcsec_pixel,1)
        ax.set_yticklabels(y_labels,size=10)
        #plt.setp(ax, xlim=custom_xlim, ylim=custom_ylim)
        ax.set_xlabel(r'arcsec',fontsize=11)
        ax.set_ylabel(r'arcsec',fontsize=11)
    return o2_err    



####################################################
# Functions to plot kinematics maps

# first component
def plot_kin(idx, f_temp, err_f_temp, v50_temp, vsig_temp, v02_temp, v98_temp, spx_x=None, spx_y=None, width=None, height=None, scale=None, spx=False, rec= False):
    figure(figsize=(5,5))
    v50_line = get_variable(idx, f_temp, err_f_temp, v50_temp)
    vsig_line = get_variable(idx, f_temp, err_f_temp, vsig_temp)
    v02_line = get_variable(idx, f_temp, err_f_temp, v02_temp)
    v98_line = get_variable(idx, f_temp, err_f_temp, v98_temp)
       
    subplot(221)
    sub_plot(v50_line, v50_min_c1, v50_max_c1,'RdBu_r', spx_x, spx_y, width, height, spx, rec,'$v_{50}\ \mathrm{(km\ s^{-1})}$', scale, lab_x=True)

    subplot(222)
    sub_plot(vsig_line, vsig_min_c1, vsig_max_c1,'jet',spx_x, spx_y, width, height, spx, rec,'$\sigma\ \mathrm{(km\ s^{-1})}$', scale, lab_x=True, lab_y=True)

    subplot(223)
    sub_plot(v02_line, v02_min_c1, v02_max_c1,'Reds',spx_x, spx_y, width, height, spx, rec,'$v_{02}\ \mathrm{(km\ s^{-1})}$', scale)

    subplot(224)
    sub_plot(v98_line, v98_min_c1, v98_max_c1,'Blues_r', spx_x, spx_y,width, height, spx, rec,'$v_{98}\ \mathrm{(km\ s^{-1})}$', scale, lab_y=True)

    plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.5,hspace=0.1)

# second component    
def plot_kin_c2(idx, f_temp, err_f_temp, v50_temp, vsig_temp, v02_temp, v98_temp, spx_x=None, spx_y=None, width=None, height=None, scale=None, spx=False, rec= False):
    figure(figsize=(5,5))
    v50_line = get_variable(idx, f_temp, err_f_temp, v50_temp)
    vsig_line = get_variable(idx, f_temp, err_f_temp, vsig_temp)
    v02_line = get_variable(idx, f_temp, err_f_temp, v02_temp)
    v98_line = get_variable(idx, f_temp, err_f_temp, v98_temp)
       
    subplot(221)
    sub_plot(v50_line, v50_min_c2, v50_max_c2,'RdBu_r', spx_x, spx_y, width, height, spx, rec,'$v_{50}\ \mathrm{(km\ s^{-1})}$', scale, lab_x=True)

    subplot(222)
    sub_plot(vsig_line, vsig_min_c2, vsig_max_c2,'jet',spx_x, spx_y, width, height, spx, rec,'$\sigma\ \mathrm{(km\ s^{-1})}$', scale, lab_x=True, lab_y=True)

    subplot(223)
    sub_plot(v02_line, v02_min_c2, v02_max_c2,'Reds',spx_x, spx_y, width, height, spx, rec,'$v_{02}\ \mathrm{(km\ s^{-1})}$', scale)

    subplot(224)
    sub_plot(v98_line, v98_min_c2, v98_max_c2,'Blues_r', spx_x, spx_y,width, height, spx, rec,'$v_{98}\ \mathrm{(km\ s^{-1})}$', scale, lab_y=True)

    plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.5,hspace=0.1)

# third component      
def plot_kin_c3(idx, f_temp, err_f_temp, v50_temp, vsig_temp, v02_temp, v98_temp, spx_x=None, spx_y=None, width=None, height=None, scale=None, spx=False, rec= False):
    figure(figsize=(5,5))
    v50_line = get_variable(idx, f_temp, err_f_temp, v50_temp)
    vsig_line = get_variable(idx, f_temp, err_f_temp, vsig_temp)
    v02_line = get_variable(idx, f_temp, err_f_temp, v02_temp)
    v98_line = get_variable(idx, f_temp, err_f_temp, v98_temp)
       
    subplot(221)
    sub_plot(v50_line, v50_min_c3, v50_max_c3,'RdBu_r', spx_x, spx_y, width, height, spx, rec,'$v_{50}\ \mathrm{(km\ s^{-1})}$', scale, lab_x=True)

    subplot(222)
    sub_plot(vsig_line, vsig_min_c3, vsig_max_c3,'jet',spx_x, spx_y, width, height, spx, rec,'$\sigma\ \mathrm{(km\ s^{-1})}$', scale, lab_x=True, lab_y=True)

    subplot(223)
    sub_plot(v02_line, v02_min_c3, v02_max_c3,'Reds',spx_x, spx_y, width, height, spx, rec,'$v_{02}\ \mathrm{(km\ s^{-1})}$', scale)

    subplot(224)
    sub_plot(v98_line, v98_min_c3, v98_max_c3,'Blues_r', spx_x, spx_y,width, height, spx, rec,'$v_{98}\ \mathrm{(km\ s^{-1})}$', scale, lab_y=True)

    plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.5,hspace=0.1)
    
def sub_plot(var,vmin,vmax,color_map,spx_x=None, spx_y=None, width=None, height=None,spx=False, rec= False,title=None,scale=None,  lab_x=False, lab_y=False):
    imshow(var, origin='lower', cmap=color_map, vmin=vmin, vmax=vmax, interpolation='none')
    ax = gca()
    cbar = colorbar(fraction=0.05, pad=0.04)
    cbar.ax.tick_params(labelsize=9)

    plt.setp(ax, xlim=custom_xlim, ylim=custom_ylim)

    if spx == True:
        plt.plot(spx_x, spx_y, marker="x", markersize=5, markeredgecolor="magenta", markerfacecolor="magenta")
    if rec == True:
        rect = patches.Rectangle((spx_x, spx_y), width, height, linewidth=1, edgecolor='magenta', facecolor='none')
        ax.add_patch(rect) 
    ax.title.set_text(title)  
    
    ax.set_xticks([16, 36, 56])
    ax.set_yticks([30, 50, 70]) 
    
    if scale == 'kpc':
        ax.set_xticklabels(np.round(np.array([-20, 0 , 20])*pixel_scale,1),size=9)
        ax.set_yticklabels(np.round(np.array([-20, 0 , 20])*pixel_scale,1),size=9)
        label = r'kpc'
    elif scale == 'arcsec':
        ax.set_xticklabels(np.round(np.array([-20, 0 , 20])*arcsec_pixel,1),size=9)
        ax.set_yticklabels(np.round(np.array([-20, 0 , 20])*arcsec_pixel,1),size=9)
        label = r'arcsec'
    else:
        label=r'spaxel'

    if lab_x == False: 
        ax.set_xlabel(label,fontsize=11)
    if lab_y == False:  
        ax.set_ylabel(label,fontsize=11)

        
####################################################
# Function to get data of one spaxel

def get_values(idx,  f_temp, err_f_temp, v50_temp, vsig_temp, v02_temp, v98_temp,  spx_x, spx_y):
    v50_line = get_variable(idx, f_temp, err_f_temp, v50_temp)
    vsig_line = get_variable(idx, f_temp, err_f_temp, vsig_temp)
    v02_line = get_variable(idx, f_temp, err_f_temp, v02_temp)
    v98_line = get_variable(idx, f_temp, err_f_temp, v98_temp)
       
    o2 = f_temp[ion][1][spx_y][spx_x]
    o2_err = err_f_temp[ion][1][spx_y][spx_x]

    ## nat added
    # print('flux', round(f_temp[spx_y][spx_x],2))

    print('v50', round(v50_line[spx_y][spx_x],2), 'vsig', round(vsig_line[spx_y][spx_x],2), 'v02', round(v02_line[spx_y][spx_x],2), 'v98', round(v98_line[spx_y][spx_x],2))
    print('o2 flux', round(o2,4), 'o2_err', round(o2_err,4), 'SNR', round(o2/o2_err,2))
    
####################################################
# Functions to plot fits

import matplotlib.ticker as ticker

class Spaxel:
    def __init__(self, infile):
        self.xdr = readsav(infile)
        self.xdr_arr = self.xdr['struct']
        self.wave = self.xdr_arr['wave'][0]
        self.spectot = self.xdr_arr['spec'][0]
        self.specstars = self.xdr_arr['cont_dat'][0]
        self.speclines = self.xdr_arr['emlin_dat'][0]
        self.modstars = self.xdr_arr['cont_fit'][0]
        self.modlines = self.xdr_arr['emlin_fit'][0]
        self.specerr = self.xdr_arr['spec_err'][0]
        self.modtot = self.modstars + self.modlines
        if self.xdr_arr['param'] !=0:
            self.ppoff = self.xdr_arr['param'][0][0]
            self.ncomp = self.xdr_arr['param'][0][1].astype(int)
            self.specres = self.xdr_arr['param'][0][2]
            self.rem_lis = []
        else:
            self.ncomp  = 0
        print(self.ncomp)
        
    def cmplin(self, line, comp):
        if comp == 0:
            return 0
        c = 299792.458
        indices = (self.xdr_arr['parinfo'][0]['line'] == line) & (self.xdr_arr['parinfo'][0]['comp'] == comp)
        if indices[0] != -1:
            gausspar = self.xdr_arr['param'][0][indices]
            gaussparRound = np.around(gausspar, decimals = 7)
            gaussparStr = list(map(str, gaussparRound))
            units = ['', ' Ang', ' km/s']
            gaussparFinal = []
            for i in range(len(gaussparStr)):
                gaussparFinal.append(gaussparStr[i] + units[i])
            #print(gaussparFinal)
            
            if len(gausspar) > 0 and not line in self.rem_lis:
                gausspar[2] = sqrt((gausspar[2]*gausspar[1]/c)**2. + self.specres**2.)
                if gausspar[2] == 0 or gausspar[0] == 0:
                    flux = zeros(len(self.wave)) * nan
                else:
                     flux = gaussian(self.wave,*gausspar)
            else:
                flux = zeros(len(self.wave)) * nan
        return flux
    
    def plot_lines(self, ax, lines, xmin, xmax):
        ax.plot(self.wave, self.modstars + self.speclines, 'k', ds='steps-mid', lw=1.)
        col = ['','0','1','2']
        for l in lines:
            for i in range(1,self.ncomp+1):
                style = 'C' + col[i] + '--'
                ax.plot(self.wave,self.modstars + self.cmplin(l, i), style)
        ax.plot(self.wave, self.modstars + self.modlines, 'C3')
        ax.plot(self.wave, self.specerr, 'C6:', ds='steps-mid', lw=1.)
        self.limits(self.wave, xmin, xmax, 2.8)

    def plot_line_vel(self, ax, line, wave, z, vmin, vmax, lw=None):
        c = 299792.458
        wobs = wave * (1. + z)
        vel = c * (self.wave - wobs) / wobs
        ax.plot(vel, self.speclines, 'k', drawstyle='steps-mid', lw=1.)
        col = ['','0','1','2']
        for l in line:
            for i in range(1,self.ncomp+1):
                style = 'C' + col[i] + '--'
                ax.plot(vel, self.cmplin(l, i), style, lw=1.)
        ax.plot(vel, self.modlines, 'C3', lw=lw)
        self.limits(vel, vmin, vmax, 1.8)

        ax.axvline(0,color='k', ls='--', lw=0.5)
        ax.set_yticks(())
        ax.set_xlabel('(km s$^{-1}$)', fontsize=6, labelpad=-1)
        ax.set_ylabel('flux')

    def plot_lines_pdf(self, lines, xmin, xmax, str_x, str_y, pr=False):
        plt.plot(self.wave, self.modstars + self.speclines, 'k', ds='steps-mid', lw=1.)
        col = ['','0','1','2']
        for l in lines:
            for i in range(1,self.ncomp+1):
                style = 'C' + col[i] + '--'
                plt.plot(self.wave,self.modstars + self.cmplin(l, i), style)
        plt.plot(self.wave, self.modstars + self.modlines, 'C3')
        plt.plot(self.wave, self.specerr, 'C6:', ds='steps-mid', lw=1.)
        ymax = self.limits(self.wave, xmin, xmax, 1.3)
        
        if pr == True:
            plt.text(plot_limits[0]+15, ymax*0.9,  str_x + '_' + str_y, fontsize=9)
            
    def limits(self, var, _min, _max, factor):
        yi = self.speclines[(var > _min) & (var < _max)] + self.modstars[(var > _min) & (var < _max)]
        erri = self.specerr[(var > _min) & (var < _max)]
        ymax = max(factor * amax(yi), 3. * amax(erri))
        ymin =  -ymax*0.2
        plt.xlim(_min, _max)
        plt.ylim(ymin, ymax)
        return ymax
        
    def plot_line_vel_pdf(self, line, wave_o, z, vmin, vmax,       str_x, str_y, sn1=None, sn2=None, lw=None, sn = False):
        c = 299792.458
        wobs = wave_o * (1. + z)
        vel = c * (self.wave - wobs) / wobs
        plt.plot(vel, self.speclines, 'k', drawstyle='steps-mid', lw=1.)
        if self.ncomp!=0:
            col = ['','0','1','2']
            for l in line:
                for i in range(1,self.ncomp+1):
                    style = 'C' + col[i] + '--'
                    plt.plot(vel, self.cmplin(l, i), style, lw=1.)
            plt.plot(vel, self.modlines, 'C3', lw=lw)
        ymax = self.limits(vel, vmin, vmax, 1.2)
        plt.axvline(0,color='k', ls='--', lw=0.5)
        if sn == True and self.ncomp!=0:
            plt.text(100, ymax*0.9,  'SNR1='+str(round(sn1,1)), fontsize=7)
            plt.text(100, ymax*0.8,  'SNR2='+str(round(sn2,1)), fontsize=7)
        ##nat added
        plt.text(plot_limits[0]+15, ymax*0.9,  str_x + '_' + str_y, fontsize=9)
            
        
def gaussian(x, A, u, sig):
    return A * exp(-(x-u)**2. / 2. / sig**2.)

def plot_spaxel(ax, spax, lxticks=False, ltitle=False, rem=False):
    if ion_tag == 'o2': 
        spax.plot_lines(ax, ion_line, plot_limits[0], plot_limits[1])
        ax.axvline((z+1)*ion_wave_0, color='k', ls='--', lw=0.5)
        ax.axvline((z+1)*ion_wave_1, color='k', ls='--', lw=0.5)
        if ltitle:
            ax.title.set_text('[OII]3726,9')
    elif ion_tag == 'o3':
        spax.plot_lines(ax, ion_line, plot_limits[0], plot_limits[1])
        ax.axvline((z+1)*ion_wave_0, color='k', ls='--', lw=0.5)
        ax.axvline((z+1)*ion_wave_1, color='k', ls='--', lw=0.5)
        if ltitle:
            ax.title.set_text('[OIII]5007')

def plot_spaxel_pdf(spax, str_x, str_y, lxticks=False, ltitle=False, rem=False, pr= False):
        spax.plot_lines_pdf(ion_line, plot_limits[0], plot_limits[1], str_x, str_y, pr)
        plt.axvline((z+1)*ion_wave_0, color='k', ls='--', lw=0.5)
        plt.axvline((z+1)*ion_wave_1, color='k', ls='--', lw=0.5)

####################################################

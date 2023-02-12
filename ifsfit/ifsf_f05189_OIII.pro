; docformat = 'rst'
;
;+
;
; This function initializes the fitting parameters for PG1700
; 
; :Categories:
;    IFSF
;
; :Returns:
;    A structure with tags specified in INITTAGS.txt.
;
; :Params:
; 
; :Keywords:
;    initmaps: out, optional, type=structure
;      Parameters for map making.
;    initnad: out, optional, type=structure
;      Parameters for NaD fitting.
; 
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2016jan23, DSNR, created
;    
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
function ifsf_f05189_OIII,initmaps=initmaps,initnad=initnad

  bad=1d99

  gal = 'F05189'
  ncols = 74
  nrows = 97
  centcol = 37
  centrow = 49
  fitrange = [3800,5300]

  platescale = 0.29d
  bin = 1d
  outstr = 'rb'+string(bin,format='(I0)')

  iter = '4'
  dir = 'r1_OIII'
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Required pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Input file
;  infile='/Users/drupke/ifs/kcwi/cubes/'+gal+'/pg1700.fits'
;  infile='/home/gene/mrk273/NONAS/mrk273_mosaic.fits'
  infile='/raid/KCWI/AGN_outflows/mosaics/F05189_mosaic.fits'
  if ~ file_test(infile) then message,'Data cube not found.'

; Lines to fit.

;  lines = ['[NeV]3426', '[OII]3726','[OII]3729', '[NeIII]3869', 'H8',$
;          '[NeIII]3967', '[SII]4068', '[SII]4076', 'Hdelta',$
;          'Hgamma', '[OIII]4363', 'HeII4686', 'Hbeta',$
;           '[OIII]4959', '[OIII]5007', '[NI]5200']

;;;;; NAT CHANGED THIS
    lines = ['[OIII]4959', '[OIII]5007']
  nlines = n_elements(lines)

; Max no. of components.
  maxncomp = 2

; Initialize line ties, n_comps, z_inits, and sig_inits.
  ncomp = hash(lines)
  linetie = hash(lines)
  zinit_gas = hash(lines)
  siginit_gas = hash(lines)
  zinit_stars=dblarr(ncols,nrows) + 0.042884d
  foreach i,lines do begin
     linetie[i] = '[OIII]5007'
     ncomp[i] = dblarr(ncols,nrows)+maxncomp
     siginit_gas[i] = dblarr(ncols,nrows,maxncomp) + 100d
     zinit_gas[i] = dblarr(ncols,nrows,maxncomp) + 0.04283d
     siginit_gas[i,*,*,1] = 100d
     zinit_gas[i,*,*,1] = 0.042884d
  endforeach

  zsys_gas = 0.043118d ;added 0.000233d (70km/s) ;0.043385d (150km/s)  ;0.0430d ;0.042885d

;  vormap = mrdfits(infile,2,/silent)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Optional pars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Tweaked regions are around [OII], [NeIII], [OIII]

; Number of wavelength regions to re-fit
;  nregions = 1
;  ; Parameters for continuum fit
;  tweakcntfit = dblarr(ncols,nrows,3,nregions)
;  ; Default fitting order
;  tweakcntfit[*,*,2,*] = 2
; Lower wavelength for re-fit
;  tweakcntfit[*,*,0,0:nregions-1] = $
;     rebin(reform([4600],1,1,1,nregions),$
;           ncols,nrows,1,nregions)
; Upper wavelength for re-fit
;  tweakcntfit[*,*,1,0:nregions-1] = $
;     rebin(reform([5500],1,1,1,nregions),$
;           ncols,nrows,1,nregions)
; Order for re-fit
;  tweakcntfit[*,*,2,0:nregions-1] = $
;     rebin(reform([6],1,1,1,nregions),$
;           ncols,nrows,1,nregions)
;
; Parameters for emission line plotting
  linoth = strarr(1,1)   ; removed all other lines, changed (1,4) to (1,1)
  linoth[0,0] = '[OIII]5007'
  
;  argspltlin1 = {nx: 3, ny: 3,$
;                 label: ['[NeV]3426','[OII]3726,3729',$
;			 '[NeIII]3869','[NeIII]3967',$
;                         '[OIII]4363', 'HeII4686',$
;                         'Hbeta','[OIII]5007'],$
;                 wave: [3426,3727,3869,3967,4363,4686,4861,4985],$
;                 off: [[-100,100],[-100,100],[-100,100],[-100,100],$
;			[-100,100],[-100,100],[-100,100],[-100,100]],$
;                 linoth: linoth}
;;;;; changed to correspond with only 1 line 
    argspltlin1 = {nx: 1, ny: 1,$
                 label: '[OIII]5007',$
                 wave: [5007],$
                 off: [-100,100],$
                 linoth: linoth}

; Velocity dispersion limits and fixed values
  siglim_gas = [50d,1000d]
  lratfix = hash()
; 1 corresponds to n ~ 400 cm^-3; Pradhan et al. 2006, MNRAS, 366, L6
;  lratfix['[OII]3729/3726']=[1.2d,1.2d]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  init = {$
; Required pars
         fcninitpar: 'ifsf_kcwi',$
         fitran: fitrange,$
         fluxunits: 1d-16, $ ; erg/s/cm^2?
         infile: infile,$
         label: gal,$
         lines: lines,$
         linetie: linetie,$
         maxncomp: maxncomp,$
         name: 'F05189',$
         ncomp: ncomp,$
;         mapdir: '/home/gene/mrk273/NONAS/ifsfit_full/maps/',$
;         outdir: '/home/gene/mrk273/NONAS/ifsfit_full/products/',$
         mapdir: '/raid/KCWI/AGN_outflows/F05189/iter'+iter+'/',$
         outdir: '/raid/KCWI/AGN_outflows/F05189/iter'+iter+'/',$
   
         platescale: platescale,$
;         positionangle: -45d,$
         positionangle: 0d,$ 
;         minoraxispa: 127d,$ ; based on major-axis PA of -153deg from Veilleux+09
         zinit_stars: zinit_stars,$
         zinit_gas: zinit_gas,$
         zsys_gas: zsys_gas,$
; Optional pars
         datext: -1,$
         varext: 1,$
         dqext: -1,$
;         vormap: vormap,$
         argscheckcomp: {sigcut: 2d},$
         argscontfit: {siginit_stars: 50d,$
                       uselog: 1b},$
;                       refit: 1b},$
         argsinitpar: {siglim: siglim_gas},$
  ;                     lratfix: lratfix},$
;         startempfile: '/home/gene/LZIFU/LZIFU-1.1/stellar_models/'+$
;                       'gonzalezdelgado/SSPGeneva_z020.sav',$
;                       gal+'hosttemplate.xdr', $
        startempfile: '/raid/KCWI/AGN_outflows/F05189/r1_OIII/F05189_starlighttemplate.xdr', $
         argspltlin1: argspltlin1,$
;         decompose_qso_fit: 1b,$
         decompose_ppxf_fit: 1b,$
         fcncheckcomp: 'ifsf_checkcomp',$
         fcncontfit: 'ppxf',$
;         fcncontfit: 'ifsf_fitqsohost',$
         maskwidths_def: 500d,$
;         tweakcntfit: tweakcntfit,$
;         emlsigcut: 2d,$
         logfile: '/raid/KCWI/AGN_outflows/F05189/'+dir+'/'+$
                  gal+'_fitlog',$
;         batchfile: '/home/gene/IFSFIT/ifsfit-master/common/ifsf_fitloop.pro',$
;         batchdir: '/home/gene/src/idl/batch/',$
         batchfile: '/data/home/wning/ifsfit/common/ifsf_fitloop.pro',$
         batchdir: '/data/raid/KCWI/AGN_outflows/F05189/'+dir+'/',$
         siglim_gas: siglim_gas,$
         siginit_gas: siginit_gas,$
         siginit_stars: 50d, $
;         nocvdf: 1, $
         cvdf_vlimits: [-2.5e3,2.5e3], $ ; in km/s
         cvdf_vstep: 1d, $ ; in km/s
;         cutrange: [[hb_maskctran],[hg_maskctran]], $
         host: {dat_fits: '/data/raid/KCWI/AGN_outflows/F05189/'+dir+'/'+'starlight_iter1.fits'} $  
        }

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Arguments for maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
   if keyword_set(initmaps) then begin

;      contourlevels = hash()
;       contourlevels['[OIII]5007_v%50c1'] = $
;         [-200,-150,-100,-50,0,50,100,150,200]

      initmaps = {$
                  emlplot: {ftags: ['ftot','fc1','fc2'],$
                            vtags: ['vsigc1','vsigc2','v%50c1','v%50c2'],$
 ;                           radtags: ['ftot'],$
                            ftitles: ['F$\uptot$',$
                                      'F$\upc1$',$
                                      'F$\upc2$'],$
                            vtitles: ['$\sigma$$\upc1$',$
                                      '$\sigma$$\upc2$',$
                                      'v$\down50$$\upc1$',$
                                      'v$\down50$$\upc2$']},$
;                            radtitles: ['F$\uptot$']},$
                  aspectrat: (double(nrows)/double(ncols)),$
                  center_axes: [centcol,centrow],$
                  center_nuclei: [centcol,centrow],$
;                  ctradprof_psffwhm: 0.6d,$
;                  contourlevels: contourlevels,$
                  ct: {sumrange: fitrange,$
;                       sumrange_hstcomp: [5250,5525],$
                       scllim: [0,1],$
                       scllim_rad: [-3,0],$
                       stretch: 1},$
;                       fitifspeak: 1b,$
;                       fitifspeakwin_kpc: 6d},$
;                  hst: {refcoords: [722,744],$
;                        subim_big: 7d},$
;                  hstrd: {file: '/Users/drupke/ifs/gmos/ancillary/hst/'+$
;                                gal+'/pg1700_acshrc_550m.fits',$
;                          label: 'ACS-HRC/F550M',$
;                          scllim: [0.001,5],$
;                          sclargs_big: {beta: 0.05,stretch: 5},$
;                          sclargs_fov: {beta: 0.05,stretch: 5},$
;                          platescale: 0.025d,$
;                          nucoffset: [0d,0d]},$
                  rangefile: '/raid/KCWI/AGN_outflows/F05189/'+$
                             'rangefile.txt',$
                  fluxunits: 1d-16 $
                 }
   endif

   initnad = {}
                  
   return,init

end

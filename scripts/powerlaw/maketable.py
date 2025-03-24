#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from astropy.io import fits
import os

variances=np.linspace(0.01,1,10)
correlations=np.linspace(0,2,10)
lengths=[10,10]
names=["Var","Pivot"]
nmax=max(lengths)
Nvars=len(lengths)
outfil="fvar_pow.fits"
if os.path.exists(outfil):
    os.remove(outfil)


# read energies from qdp file:
energies=np.loadtxt("dump_baseline.txt",skiprows=3)[:,0]
energies_err=np.loadtxt("dump_baseline.txt",skiprows=3)[:,1]

i=0

spectra=[]
parray=[]

# ax=pl.subplot(111)
# ax.set_xscale("log")
# ax.set_yscale("log")

for var in variances:
    for cor in correlations:
        filename="grid_spectra/plsim_var%s_cor%s.txt" % (str(var),str(cor))
        spectra.append(np.loadtxt(filename))
        parray.append([var,cor])
        # pl.plot(energies,spectra[-1]/energies)

# pl.show()
spectra=np.array(spectra)
# test_spectra=np.ones(spectra.shape)
parray=np.array(parray)



prihdu = fits.PrimaryHDU()
prihd =prihdu.header
prihd.extend([('MODLNAME','FVAR_POW'),('MODLUNIT','PHOTONS/CM2/S'),\
                  ('REDSHIFT',True),('ADDMODEL',True),('HDUCLASS','OGIP'),\
                  ('HDUCLAS1','XSPEC TABLE MODEL'),('HDUVERS','1.0.0'),\
                  ('AUTHOR','M L PARKER'),\
                  ('COMMENT','')])

pcnames = ['NAME','METHOD','INITIAL','DELTA','MINIMUM','BOTTOM',\
               'TOP','MAXIMUM','NUMBVALS','VALUE']
pcformats = ['12A','J','E','E','E','E','E','E','J','%sE' % nmax]

p1=["Var",0,0.5,0.1,0.01,0.01,1,1,nmax,variances]
p2=["Corr",0,0.5,0.1,0,0,2,2,nmax,correlations]


pars=[p1,p2]

parcols=[]
for c in range(0,len(pars[0])):
    col=[]
    for p in range(0,Nvars):
        par = pars[p]
        col.append(par[c])
    # print(col)
    parcols.append(fits.Column(name=pcnames[c],format=pcformats[c],array=col))


pcdefs = fits.ColDefs(parcols)
partb = fits.BinTableHDU.from_columns(pcdefs)
partb.name='Parameters'
parhd = partb.header
parhd.extend([('NINTPARM',Nvars),('NADDPARM',0),('HDUCLASS','OGIP'),\
                  ('HDUCLAS1','XSPEC TABLE MODEL'),\
                  ('HDUCLAS2','PARAMETERS'),('HDUVERS','1.0.0')])


energ_lo = fits.Column(name='ENERG_LO', format='E', array=energies-energies_err)
energ_hi = fits.Column(name='ENERG_HI', format='E', array=energies+energies_err)
binwidth=energies_err*2.

energtb = fits.BinTableHDU.from_columns([energ_lo,energ_hi])
energtb.name = 'Energies'
energhd = energtb.header
energhd.extend([('HDUCLASS','OGIP'),('HDUCLAS1','XSPEC TABLE MODEL'),\
                  ('HDUCLAS2','ENERGIES'),('HDUVERS','1.0.0')])




parcol = fits.Column(name = 'PARAMVAL',format='%sE' %Nvars ,array = parray)
speccol = fits.Column(name = 'INTPSPEC',format='%sE' % len(spectra[0]),\
                        unit='photons/cm2/s/keV',array = spectra*binwidth)

spectb = fits.BinTableHDU.from_columns([parcol,speccol])
spectb.name = 'Spectra'
spechd = spectb.header
spechd.extend([('HDUCLASS','OGIP'),('HDUCLAS1','XSPEC TABLE MODEL'),\
                  ('HDUCLAS2','MODEL SPECTRA'),('HDUVERS','1.0.0')])


thdulist = fits.HDUList([prihdu, partb, energtb, spectb])

thdulist.writeto(outfil)

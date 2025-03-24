#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from astropy.io import fits
import os

# variances=np.linspace(0.01,1,10)
# correlations=np.linspace(0,2,10)
correlations=np.linspace(0,1,10)
xis=np.linspace(3,5,3)
nhs=np.logspace(-1,1,5)
lengths=[10,5,5]
nmax=max(lengths)
Nvars=len(lengths)
outfil="fvar_ufo.fits"
if os.path.exists(outfil):
    os.remove(outfil)


# read energies from qdp file:
energies=np.loadtxt("dump_baseline.txt",skiprows=3)[:,0]
energies_err=np.loadtxt("dump_baseline.txt",skiprows=3)[:,1]

i=0

spectra=[]
parray=[]

ax=pl.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")

for xi in xis:
    for nh in nhs:
        for cor in correlations:

            filename="grid_spectra/ufosim_xi%s_nh%s_cor%s.txt" % (str(xi),str(nh),str(cor))
            spectra.append(np.loadtxt(filename))
            parray.append([xi,nh,cor])
            pl.plot(energies,spectra[-1])


pl.show()
ref_spec=np.loadtxt("../powerlaw/grid_spectra/plsim_var0.2_cor0.3.txt")

spectra=np.array(spectra)
specmax=np.max(spectra,axis=0)
specmin=np.min(spectra,axis=0)
specmed=np.median(spectra,axis=0)
# print(specmed)

# ax2=pl.subplot(111)
# ax2.set_xscale("log")
# pl.plot(energies,specmed)
# pl.plot(energies,ref_spec)
# pl.show()
spectra=spectra/ref_spec
# test_spectra=np.ones(spectra.shape)
parray=np.array(parray)
# exit()



prihdu = fits.PrimaryHDU()
prihd =prihdu.header
prihd.extend([('MODLNAME','FVAR_UFO'),('MODLUNIT','PHOTONS/CM2/S'),\
                  ('REDSHIFT',True),('ADDMODEL',False),('HDUCLASS','OGIP'),\
                  ('HDUCLAS1','XSPEC TABLE MODEL'),('HDUVERS','1.0.0'),\
                  ('AUTHOR','M L PARKER'),\
                  ('COMMENT','')])

pcnames = ['NAME','METHOD','INITIAL','DELTA','MINIMUM','BOTTOM',\
               'TOP','MAXIMUM','NUMBVALS','VALUE']
pcformats = ['12A','J','E','E','E','E','E','E','J','%sE' % nmax]


xi_array=np.zeros(nmax)
for i,xi in enumerate(xis):
    xi_array[i]=xi

nh_array=np.zeros(nmax)
for i,nh in enumerate(nhs):
    nh_array[i]=nh
print(xi_array)

p1=["logxi",0,3.5,0.1,2,2,5,5,3,xi_array]
p2=["nH",1,1,0.1,0.1,0.1,10,10,5,nh_array]
p3=["cor",0,0.5,0.1,0,0,1,1,nmax,correlations]


pars=[p1,p2,p3]

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
                        unit='photons/cm2/s/keV',array = spectra)

spectb = fits.BinTableHDU.from_columns([parcol,speccol])
spectb.name = 'Spectra'
spechd = spectb.header
spechd.extend([('HDUCLASS','OGIP'),('HDUCLAS1','XSPEC TABLE MODEL'),\
                  ('HDUCLAS2','MODEL SPECTRA'),('HDUVERS','1.0.0')])


thdulist = fits.HDUList([prihdu, partb, energtb, spectb])

thdulist.writeto(outfil)

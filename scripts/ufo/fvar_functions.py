
import numpy as np
import os
import matplotlib.pyplot as pl
from xspec import *
from scipy.ndimage.filters import gaussian_filter1d as gsmooth
from matplotlib.ticker import *
from scipy.interpolate import UnivariateSpline as USpline


os.environ['RELLINE_TABLES']='/Users/mlparker/programs/xspec/relxill/'
AllModels.lmod('relxill','/Users/mlparker/programs/xspec/relxill/')
AllModels.lmod('absmodel','/Users/mlparker/programs/xspec/tbnew/')

def parse_parfile(filename):
    infile=open(filename,"r")
    print("Reading parameters to vary:")
    pars=[]
    for row in infile:
        if row[0]!="!":
            splitrow=row.split()
            pars.append((int(splitrow[0]),float(splitrow[1])))
            print("Parameter %s will be varied with a sigma of %s" % (splitrow[0],splitrow[1]))
    if len(pars)==0:
        print("WARNING: No parameters read in parse_parfile from file",filename)
    return pars

def read_xcm(xcmfilename):
    '''Read data, model etc. from an xspec .xcm file
    Note that this is not complete. Check the list below.
    '''

    cwd = os.getcwd()

    print('\nReading .xcm file...')
    xcmfile = open(xcmfilename)

    model_index = -1

    pars = {}
    newpars = {}
    n_datasets=1
    for i, row in enumerate(xcmfile):
        temp = row.strip().split()

        # print row
        if len(temp) > 0:
            # Change directory
            if temp[0] == 'cd':
                os.chdir(temp[1])

            # Set abundances and cross-sections
            elif temp[0] == 'abund':
                Xset.abund = temp[1]
            elif temp[0] == 'xsect':
                Xset.xsect = temp[1]

            # Load data
            elif temp[0] == 'data':
                print('Loading data:', ' '.join(temp[1:]))
                AllData(' '.join(temp[1:]))
                n_datasets = int(temp[1][0])

            # Ignore specified channels
            elif temp[0] == 'ignore':
                AllData.ignore(' '.join(temp[1:]))

            # Find model definition
            elif temp[0] == 'model':
                print('Loading model:', ' '.join(temp[1:]))
                model_index = i
                modelstr = row.strip()[7:]
                model = Model(modelstr)

            elif model_index > 0 and i > model_index:
                if temp[0] == 'newpar':
                    parnum = int(temp[1])
                    newpars[parnum] = ' '.join(temp[2:])
                else:
                    pars[i - model_index] = ' '.join(temp)

    n_pars = len(pars)
    if n_pars % n_datasets != 0:
        print(n_pars, n_datasets, n_pars % n_datasets)
        print("Number of parameters should be a multiple of the number of datasets. This isn't.")

    else:
        n_pars_per_dset = n_pars / n_datasets
        for p in pars:
            # print(p, pars[p])
            mnum = int(1 + (p - 1) / n_pars_per_dset)
            pnum = int(p - (mnum - 1) * n_pars_per_dset)
            # print(mnum, pnum)
            AllModels(mnum).setPars({pnum: pars[p]})
            # print(pnum, pars[p])
        for p in newpars:
            # print p, newpars[p]
            mnum = int(1 + (p - 1) / n_pars_per_dset)
            pnum = int(p - (mnum - 1) * n_pars_per_dset)

            AllModels(mnum).setPars({pnum: newpars[p]})

    return pars, n_datasets


def vary_parameters(variable_pars,N, cor, xi, nh):
    print("Varying parameters and calculating spectra...")
    Plot.device = "/xw"
    Plot("eemodel")

    gamma_cor=0.3
    gamma_offset=2+12.*gamma_cor
    gamma_par=7
    pl_flux_par=6
    xi_par=1
    xi_offset=xi+12*cor
    nh_par=2

    initial_pars=[]
    for par in variable_pars:
        initial_pars.append(AllModels(1)(par[0]).values)


    spectra=[]

    # Temp fudge:
    # print(frac)
    # print(initial_pars[0])
    # print(bb_flux)
    AllModels(1).setPars({xi_par:xi})
    AllModels(1).setPars({nh_par:nh})
    for i in range(0,N):
        for j,par in enumerate(variable_pars):
            pmin=initial_pars[j][2]
            pmax=initial_pars[j][5]
            pmean=initial_pars[j][0]

            newval=np.random.normal(pmean,par[1])
            newval=max(newval,pmin)
            newval=min(newval,pmax)
            # print(newval)

            AllModels(1).setPars({par[0]:newval})

        # Any extra fudges go here:
        # Gamma fudge:
        pl_flux=AllModels(1)(pl_flux_par).values[0]
        new_gamma=gamma_cor*pl_flux+gamma_offset
        if new_gamma<-3:
            AllModels(1).setPars({gamma_par:-3})
        else:
            AllModels(1).setPars({gamma_par:new_gamma})
        #Ionization fudge
        new_xi=cor*pl_flux+xi_offset
        if new_xi<2:
            AllModels(1).setPars({xi_par:2})
        elif new_xi>5:
            AllModels(1).setPars({xi_par:5})
        else:
            AllModels(1).setPars({xi_par:new_xi})


        Plot()
        spectra.append(Plot.model())

    Es=Plot.x()

    for j,par in enumerate(variable_pars):
        AllModels(1).setPars({par[0]:initial_pars[j][0]})

    return np.array(spectra),Es

def plot_spectra(spectra,energies):
    spec_fig=pl.figure()
    ax=pl.subplot(111)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.4,10)
    for spectrum in spectra:
        pl.plot(energies,gsmooth(spectrum,20),color="k",lw=1,alpha=0.1)

    pl.plot(energies,gsmooth(np.mean(spectra,axis=0),20),color='r')

    pl.savefig("spectra.pdf",bbox_inches="tight")


def calc_variance(spectra,energies):
    # print(wna_data.shape)
    var=np.var(spectra,axis=0)
    meanspec=np.mean(spectra,axis=0)
    Fvar=(var/meanspec**2)**0.5
    Fvar[meanspec<=0]=0
    return Fvar

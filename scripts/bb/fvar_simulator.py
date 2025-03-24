#!/usr/bin/env python

from fvar_functions import *

#### Settings. Change these as appropriate.
xcm_filename="baseline.xcm"

var=0.2
cor=0.3
kTs=np.logspace(-2,0,10)
fractions=np.linspace(1e-3,2,10)

# Read XCM file and set up model
read_xcm(xcm_filename)

i=0
for kT in kTs:
    for frac in fractions:

        filename="grid_spectra/bbsim_kT%s_frac%s.txt" % (str(kT),str(frac))
        i+=1

        if not os.path.exists(filename):
            print("Calculating spectra for step %s of %s" % (str(i),str(len(kTs)*len(fractions))))

            variable_pars=[(3,var)]

            # Vary the variable parameters, return spectra
            spectra,energies=vary_parameters(variable_pars,1000,cor,kT,frac)

            # plot_spectra(spectra,energies)
            Fvar=calc_variance(spectra,energies)

            np.savetxt(filename,Fvar)

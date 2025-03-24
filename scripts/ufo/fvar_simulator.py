#!/usr/bin/env python

from fvar_functions import *

#### Settings. Change these as appropriate.
xcm_filename="baseline.xcm"

var=0.2
correlations=np.linspace(0,1,10)
xis=np.linspace(3,5,3)
nhs=np.logspace(-1,1,5)

# Read XCM file and set up model
read_xcm(xcm_filename)

i=0
for xi in xis:
    for nh in nhs:
        for cor in correlations:

            filename="grid_spectra/ufosim_xi%s_nh%s_cor%s.txt" % (str(xi),str(nh),str(cor))
            i+=1

            if not os.path.exists(filename):
                print("Calculating spectra for step %s of %s" % (str(i),str(len(xis)*len(nhs)*len(correlations))))

                variable_pars=[(6,var)]

                # Vary the variable parameters, return spectra
                spectra,energies=vary_parameters(variable_pars,1000,cor,xi,nh)

                # plot_spectra(spectra,energies)
                Fvar=calc_variance(spectra,energies)

                np.savetxt(filename,Fvar)

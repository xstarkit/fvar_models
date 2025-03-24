#!/usr/bin/env python

from fvar_functions import *

#### Settings. Change these as appropriate.
xcm_filename="baseline.xcm"

var=0.2
cor=0.3
densities=np.linspace(15,19,5)
correlations=np.linspace(0,1,10)
fractions=np.linspace(0,5,10)
xis=np.linspace(1,3,3)

# Read XCM file and set up model
read_xcm(xcm_filename)

i=0
for density in densities:
    for frac in fractions:
        for cor in correlations:
            for xi in xis:
                filename="grid_spectra/refsim_dens%s_frac%s_cor%s_xi%s.txt" % (str(density),str(frac),str(cor),str(xi))

                i+=1

                if not os.path.exists(filename):
                    print("Calculating spectra for step %s of %s" % (str(i),str(len(densities)*len(fractions)*len(correlations)*len(xis))))

                    variable_pars=[(20,var)]

                    # Vary the variable parameters, return spectra
                    spectra,energies=vary_parameters(variable_pars,1000,cor,density,frac,xi)

                    # plot_spectra(spectra,energies)
                    Fvar=calc_variance(spectra,energies)

                    np.savetxt(filename,Fvar)

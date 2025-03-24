#!/usr/bin/env python

from fvar_functions import *

#### Settings. Change these as appropriate.
xcm_filename="baseline.xcm"

# variances=np.linspace(0.01,1,10)
# correlations=np.linspace(0,2,10)
variances=[0.2]
correlations=[0.3]

# Read XCM file and set up model
read_xcm(xcm_filename)

i=0
for var in variances:
    for cor in correlations:

        filename="grid_spectra/plsim_var%s_cor%s.txt" % (str(var),str(cor))
        i+=1

        if not os.path.exists(filename):
            print("Calculating spectra for step %s of %s" % (str(i),str(len(variances)*len(correlations))))

            variable_pars=[(3,var)]

            # Vary the variable parameters, return spectra
            spectra,energies=vary_parameters(variable_pars,1000,cor)

            # plot_spectra(spectra,energies)
            Fvar=calc_variance(spectra,energies)

            np.savetxt(filename,Fvar)

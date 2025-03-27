## XSPEC excess variance models

Models are calculated with the Monte-Carlo method described in [Parker et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1363P/abstract). 

These models are simple, cover a limited parameter space, and rely on basic assumptions of variability. Use at your own risk.

Important! Due to a typo, the parameter labelled var should instead be labelled sigma. This will be corrected in future models.


### Guidelines for use

 - For models with a normalisation, fix this at 1, otherwise the results are not valid.
 - No Galactic absorption is required, or any other constant multiplicative component.
 - Do not use these models with count spectra. That’s not what they’re for.
 - Because the models are relatively unsophisticated compared to conventional X-ray models, we recommend adding a small 1-2% systematic error. This is probably not necessary for low quality RMS spectra.
 - Note that some versions of xspec have a bug when adding systematic error that may cause them to not fit. This can be avoided by running an initial fit without systematic error, before adding it in and re-fitting.
 - For models derived from Spex models (for example, fvar_ufo), the unit convention follows that of Spex (SI units) rather than xspec (cgs).
 - Regular response matrices are not suitable for use with these models - variance will not be redistributed in the same way as counts. For now, we recommend approximating the instrumental resolution with a Gaussian smoothing.
 - Cite [Parker et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1363P/abstract).


### New models

New models will be added to this directory as we generate them. If you need a specific model that we do not provide, please get in touch.

### Contact

mparker@sciops.esa.int (expired!)


### fvar_pow.fits ([Parker et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1363P/abstract))

Basic power-law variance model. Additive model, which other models can be applied to.
Note that the normalization should be fixed at 1, otherwise the variance parameter is meaningless.

#### Parameters

 - `var`: The variance of log(flux) of the powerlaw, between 0.3 and 10 keV.
 - `cor`: The strength of the correlation between the powerlaw flux and the photon index.


### fvar_bbdamp.fits ([Parker et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1363P/abstract))

Damping from a constant black-body. Used to model the less variable soft excess in AGN. Multiplicative model, should be used to multiply an additive model (such as fvar_pow).

#### Parameters:
 - `kT`: Black-body temperature
 - `frac`: The fraction of the 0.3-10 keV flux attributable to the black-body (assuming the only other component is a powerlaw)


### fvar_refdamp.fits ([Parker et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1363P/abstract))
Damping from a relativistic reflection spectrum. Multiplicative model, should be used to multiply an additive model (such as fvar_pow).

#### Parameters:
 - `density`: log(density) of the accretion disk 
 - `frac`: The fraction of the 0.3-10 keV flux attributable to the reflection spectrum (assuming the only other component is a powerlaw)
 - `cor`: The strength of the correlation between the powerlaw flux and the reflection flux.


### fvar_ufo.fits ([Parker et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1363P/abstract), [Härer et al., 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.4506H/abstract))

Enhancement of the variability from absorption lines, where the ionization of the absorbing gas is driven by the X-ray continuum flux. Härer et al. introduced a version where the velocity of the gas also correlates with the continuum. This is currently only available as a legacy model, as we have not yet regenerated it with the new sampling algorithm. Multiplicative model, should be used to multiply an additive model (such as fvar_pow).

#### Parameters:
 - `logxi`: Mean log ionisation of the UFO gas.
 - `nH`: Column density of the UFO gas in units of $10^{24}$/cm^2
 - `cor`: Correlation between the UFO ionisation and the powerlaw flux


### fvar_wa.fits (Danhaive et al., in prep)

This is the same model as fvar_ufo, but covering a lower ionisation range and lower column densities, appropriate for warm absorbers. Multiplicative model, should be used to multiply an additive model (such as fvar_pow).


### fvar_pidamp.fits ([Parker et al., 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.1798P/abstract))

Damping of the variability produced by ~constant photoionised emission lines, as observed in 1H0707 and NGC 4051 (Parker et al., 2021, Danhaive et al., in prep.) Multiplicative model, should be used to multiply an additive model (such as fvar_pow).

#### Parameters:
 - `frac`: The fraction of the average 0.3-10 keV flux attributable to the photoionised emission (assuming the only other component is a powerlaw)
 - `xi`: Log ionisation of the photoionised emission.


### fvar_pcov.fits ([Parker et al., 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.1798P/abstract))

Variability from partial covering absorption independent of the continuum. Varies in covering fraction, with an optional parameter allowing for a correlation between covering fraction and column density.

#### Parameters:
 - `nh`: column density of the absorber in $10^{24}$/cm^2
 - `fcov`: mean covering fraction (note that with large var at the limits of 0 and 1 this will not be the true mean, as the distribution will be truncated at 0 and 1)
 - `var`: variance of the covering fraction
 - `nh_cor`: correlation between the column density and the covering fraction
 - `xi`: log ionisation of the absorber


### fvar_xildamp.fits ([Parker et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1363P/abstract))

Damping from neutral reflection. This mainly produces a negative feature at 6.4 keV.

#### Parameters:
 - `frac`: The fraction of the average 0.3-10 keV flux attributable to distant reflection (assuming the only other component is a powerlaw)

### Original Link

https://web.archive.org/web/20250116024332/https://www.michaelparker.space/variance-models

# LIRA-Ising: Estimating the boundaries of irregularly shaped sources

Bayesian algorithm to identify the boundaries of irregularly shaped sources. First, LIRA (pylira) needs to be run to obtain the posterior pixelwise probability distributions 
of the source intensity that properly account for known structures, astrophysical background, and the effect of the telescope point spread function.  Then this algorithm  uses Ising model to group 
pixels with similar intensities into cohesive regions corresponding to background and source. The boundary is derived on the basis of the most likely aggregation of pixels into the source region. 

# QG_DNS
A spectral code for two-equal-layer QG in a periodic domain

The matlab script Driver.m runs the simulation. 
All the parameters are set in the script Initialize.m
QG_RHS.m computes the advection terms (including beta) and bottom friction. It does not compute the hyperdiffusive PV dissipation term.
Versions of this code were used, for example, in the following two papers
I Grooms and AJ Majda, "Stochastic superparameterization in quasigeostrophic turbulence" J Comput Phys, 2014
I Grooms and L Zanna, "A note on 'Toward a stochastic parameterization of mesoscale eddies'" Ocean Modelling, 2017

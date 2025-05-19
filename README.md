# QG_DNS
A spectral code for both two-equal-layer QG and two-surface QG in a periodic domain

The matlab script Driver.m runs the simulation. 

All the parameters are set in the script Initialize.m

QG_RHS.m computes the advection terms (including beta) and bottom friction. It does not compute the hyperdiffusive PV dissipation term.

Versions of this code were used, for example, in the following papers

I Grooms and AJ Majda, "Stochastic superparameterization in quasigeostrophic turbulence" J Comput Phys, 2014

I Grooms and L Zanna, "A note on 'Toward a stochastic parameterization of mesoscale eddies'" Ocean Modelling, 2017

JB Weiss and I Grooms, "Assimilation of ocean sea-surface height observations of mesoscale eddies." Chaos: An Interdisciplinary Journal of Nonlinear Science, 2017

A Chen et al., "Comparing Eddy‐Permitting Ocean Model Parameterizations via Lagrangian Particle Statistics in a Quasigeostrophic Setting" JGR Oceans, 2018

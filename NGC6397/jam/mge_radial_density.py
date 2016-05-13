"""
#############################################################################

Copyright (C) 2014-2016, Michele Cappellari
E-mail: michele.cappellari_at_physics.ox.ac.uk

Updated versions of the software are available from my web page
http://purl.org/cappellari/software

If you have found this software useful for your research,
I would appreciate an acknowledgment to the use of the
"JAM modelling package of Cappellari (2008)"
or a reference to Note 11 of Cappellari et al. (2015).

This software is provided as is without any warranty whatsoever.
Permission to use, for non-commercial purposes is granted.
Permission to modify for personal or internal use is granted,
provided this copyright and disclaimer are included unchanged
at the beginning of the file. All other rights are reserved.

#############################################################################

NAME: 
    mge_radial_density
    
PURPOSE:
    Calculates the spherically-averaged density of an axisymmetric MGE
    model, at a given set of radii, using the equation in Note 11 of
    Cappellari et al. (2015) http://adsabs.harvard.edu/abs/2015ApJ...804L..21C
  
INPUT PARAMETERS:
    SURF: Peak surface brightness of each Gaussian in Lsun/pc**2
    SIGMA: Gaussian dispersion in arcsec
    QOBS: Observed axial ratio of each Gaussian
    INC: Inclination of the MGE model, with i=90 being edge-on
    RAD: radius for which the density is requested (scalar or vector)
    
OPTIONAL INPUT KEYWORDS:
    DISTANCE: Galaxy distance in Mpc (Default 10pc = 1e-5 Mpc)

OUTPUTS:
    DENS: spherical density calculated at the input radii RAD

MODIFICATION HISTORY:
    V1.0.0: Written and tested. Michele Cappellari, Oxford, 8 April 2014
    V1.0.1: All input distances are converted to pc. MC, Oxford, 15 May 2015
    V1.0.2: Updated documentation. MC, Oxford, 2 March 2016

"""

import numpy as np
from scipy.special import erf

def mge_radial_density(surf, sigma, qobs, inc, rad, distance=1e-5):

    if not (surf.size == sigma.size == qobs.size):
        raise ValueError("The MGE components do not match")

    rad = np.atleast_1d(rad)
    pc = distance*np.pi/0.648  # Constant factor to convert arcsec --> pc
    sigma = sigma*pc
    rad = rad*pc

    # total luminosity each Gaussian (eq.12 of Cappellari 2008)
    lum = 2.*np.pi*surf*qobs*sigma**2

    # Axisymmetric deprojection (eq.14 of Cappellari 2008, MNRAS, 390, 71)
    inc = np.radians(inc)
    qintr = qobs**2 - np.cos(inc)**2
    if np.any(qintr <= 0):
        raise ValueError('Inclination too low q < 0')
    qintr = np.sqrt(qintr)/np.sin(inc)
    qintr = qintr.clip(0.001, 0.999) # Clip to avoid division by zero in dens
    
    # Average density rho(r) on spherical surfaces, without circularizing MGE.
    # Ths is the integral of a 3-dim axisymmetric Gaussian over the surface
    # of a sphere, divided by the surface of that sphere (4pi*r^2).
    # This is the formula in footnote 11 of Cappellari et al. (2015, ApJL)
    e = np.sqrt(1. - qintr**2)
    rad = rad[:, np.newaxis]
    dens = np.sum(lum*np.exp(-0.5*(rad/sigma)**2)
                  * erf(e*rad/(np.sqrt(2)*qintr*sigma))
                  / (4*np.pi*e*rad*sigma**2), axis=1)
   
    return dens

#------------------------------------------------------------------------------ 


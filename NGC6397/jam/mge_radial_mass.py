"""
#############################################################################

Copyright (C) 2011-2016, Michele Cappellari
E-mail: michele.cappellari_at_physics.ox.ac.uk

Updated versions of the software are available from my web page
http://purl.org/cappellari/software

If you have found this software useful for your research,
I would appreciate an acknowledgment to the use of the
"JAM modelling package of Cappellari (2008)"
or an explicit reference to the equations given below.

This software is provided as is without any warranty whatsoever.
Permission to use, for non-commercial purposes is granted.
Permission to modify for personal or internal use is granted,
provided this copyright and disclaimer are included unchanged
at the beginning of the file. All other rights are reserved.

#############################################################################

MODIFICATION HISTORY:
    V1.0.0: Michele Cappellari, Oxford, 12 January 2011
    V2.0.0: Translated from IDL into Python. MC, Oxford, 21 February 2014
    V2.0.1: Cleaned up code and documentation. MC, Oxford, 23 October 2015
    V2.0.3: Check input parameters. MC, Oxford, 2 March 2016

"""

import numpy as np
from scipy.special import erf

def mge_radial_mass(surf, sigma, qobs, inc, rad, distance=1e-5):
    """
    Calculates the analytic mass or luminosity of an axisymmetric MGE within
    a sphere of a given radius. This function uses the analytic integral

      M(R) = 4*Pi*Integral[r^2*rhoTot(r), {r, 0, R}]

    where the spherical density rhoTot(r) is given in
    Note 11 of Cappellari et al. (2015, ApJL, 804, L21)
    http://adsabs.harvard.edu/abs/2015ApJ...804L..21C

    :param surf: Peak surface brightness of each Gaussian in Lsun/pc**2
    :param sigma: Gaussian dispersion in arcsec
    :param qobs: Observed axial ratio of each Gaussian
    :param inc: inclination of the MGE in degrees
    :param rad: radius of the sphere in arcsec (scalar or vector)
    :param distance: Galaxy distance in Mpc
    :return: mass of the MGE within spheres of radii rad (same size as rad)
    """

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
    qintr = qintr.clip(0.001, 0.999) # Clip to avoid division by zero in mass_r
    
    # Compute analytic integral of MGE inside a sphere
    e = np.sqrt(1. - qintr**2)
    h = rad[:, np.newaxis]/(np.sqrt(2.)*sigma*qintr)
    mass_r = np.sum(lum*(erf(h) - np.exp(-(h*qintr)**2)*erf(h*e)/e), axis=1)
   
    return mass_r

#------------------------------------------------------------------------------

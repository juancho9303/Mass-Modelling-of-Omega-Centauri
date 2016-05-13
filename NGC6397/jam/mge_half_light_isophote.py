"""
#############################################################################

Copyright (C) 2010-2016, Michele Cappellari
E-mail: cappellari_at_astro.ox.ac.uk

Updated versions of the software are available from my web page
http://purl.org/cappellari/software

If you have found this software useful for your research,
I would appreciate an acknowledgment to the use of the
"JAM modelling package of Cappellari (2008)"
or an explicit reference to the equation given below.

This software is provided as is without any warranty whatsoever.
Permission to use, for non-commercial purposes is granted.
Permission to modify for personal or internal use is granted,
provided this copyright and disclaimer are included unchanged
at the beginning of the file. All other rights are reserved.

#############################################################################

NAME:
    MGE_HALF_LIGHT_ISOPHOTE

PURPOSE:
    Computes the half-light radius, the  major axis and ellipticity of
    the MGE isophote containing 1/2 of the analytic MGE total light.

    This procedure implements the steps (i)-(iv) described above
    equation (12) in Cappellari et al. (2013, MNRAS, 432, 1709)
    http://adsabs.harvard.edu/abs/2013MNRAS.432.1709C

CALLING SEQUENCE:
    reff, re_maj, eps, lum_tot = \
        mge_half_light_isophote(surf, sigma, qObs, distance=distance, npix=npix)

INPUT PARAMETERS:
    SURF: Peak surface brightness of each Gaussian in Lsun/pc**2
    SIGMA: Gaussian dispersion in arcsec
    QOBS: Observed axial ratio of each Gaussian
    
OPTIONAL INPUT KEYWORDS:
    DISTANCE: Galaxy distance in Mpc. This is only required to obtain
        the total luminosity of the MGE model in proper units.
      - If the distance is not given, 10pc is assumed. In this case the
        following expression gives the galaxy apparent total magnitude:
        mag = sunMag - 2.5*alog10(lum_tot)
    NPIX: Number of pixels used to construct one quadrant of the MGE
        image used to find the half-light isophotes (default NPIX=500).

OUTPUTS:
    REFF: Effective (projected half-light) radius Re. This is the
        "circularized" Re=sqrt(Area/pi), where the area is the one
        of the isophote containing 1/2 of the total MGE light.
        This is in the same units as SIGMA (typically arcsec).
    RE_MAJ: Major axis (largest size) of the isophote containing 
        1/2 of the total MGE light.
    EPS: Ellipticity of the isophote containing 1/2 of the total MGE
        light, computed from the inertia tensor inside the isophote.
    LUM_TOT: The total luminosity in solar luminosities if the 
        optional distance in Mpc is provided.

MODIFICATION HISTORY:
    V1.0.0: Written and tested. Michele Cappellari, Oxford, 24 April 2010
    V1.0.1: Use major axis fluxes as reference for isophotes.
      MC, Oxford, 27 July 2011
    V2.0.0: Translated from IDL into Python. MC, Oxford, 21 February 2014
    V2.1.0: Major speedup using histogram. Updated documentation.
      MC, Oxford, 11 March 2016

"""

import numpy as np

##############################################################################

def mge_half_light_isophote(surf, sigma, qObs, distance=1e-5, npix=500):

    if not (surf.size == sigma.size == qObs.size):
        raise ValueError("The MGE components do not match")

    pc = distance*np.pi/0.648  # Constant factor to convert arcsec --> pc
    lum_tot = np.sum(2*np.pi*surf*(sigma*pc)**2*qObs)  # Total MGE luminosity
     
    # Create image from MGE model. Only compute one quadrant for symmetry
    dx = np.max(sigma)  # Re cannot be larger than this
    scale = dx/npix    # image scale in arcsec/pix
    x2 = np.linspace(scale/2, dx - scale/2, npix)**2  # open interval to use bi-symmetry
    xx2, yy2 = np.meshgrid(x2, x2)
    image = np.zeros_like(xx2)
    for su, si, qo in zip(surf, sigma, qObs):  # Construct MGE model
        image += su*np.exp(-0.5/si**2*(xx2 + yy2/qo**2))

    # Find enclosed light inside isophotes defined by the MGE flux for the
    # pixels on the major axis, then interpolate to find half-light isophote
    h = np.histogram(image, bins=image[0, ::-1], weights=image)[0]
    lum_r = np.cumsum(h[::-1])
    surf_e = np.interp(lum_tot/8/(pc*scale)**2, lum_r, image[0, 1:])

    mask = image >= surf_e  # Pixels inside half-light isophote
    reff = np.sqrt(4*mask.sum()/np.pi)*scale  # Radius of same-area circle
    re_maj = np.sqrt(np.max(xx2[mask]))   # Major axis of half-light isophote

    # Compute coefficients of the diagonal moment of inertia tensor
    # using equation (12) of Cappellari et al. (2013, MNRAS, 432, 1709)
    x2 = np.sum(image[mask]*xx2[mask])
    y2 = np.sum(image[mask]*yy2[mask])
    eps = 1. - np.sqrt(y2/x2)
    
    return reff, re_maj, eps, lum_tot

##############################################################################

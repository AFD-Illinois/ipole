import h5py
import numpy as np
import ehtim as eh
from scipy.interpolate import RegularGridInterpolator

def argb2_interp(im, minR, maxR, minphi, maxphi, narr, deg=True):
    npix = im.xdim
    fov_muas = im.fovx()/eh.RADPERUAS
    iarr = np.reshape(im.imvec,(npix,npix))
    qarr = np.reshape(im.qvec,(npix,npix))
    uarr = np.reshape(im.uvec,(npix,npix))
    pxi = (np.arange(npix))/npix-0.5
    pxi += (pxi[1]-pxi[0])/2
    pxj = np.copy(pxi)
    mui = pxi*fov_muas
    muj = pxj*fov_muas
    MUI,MUJ = np.meshgrid(mui,muj)
    Reven = np.linspace(minR,maxR,narr)
    phieven = np.arange(minphi,maxphi,2*np.pi/narr)
    
    interp_I = RegularGridInterpolator((mui, muj), iarr)
    interp_Q = RegularGridInterpolator((mui, muj), qarr)
    interp_U = RegularGridInterpolator((mui, muj), uarr)
    
    betavals = []
    unitcoords = np.array([[-np.cos(phi),-np.sin(phi)] for phi in phieven])
    for R in Reven:
        coords = R*unitcoords
        Phere = interp_Q(coords)+1j*interp_U(coords)
        integrand = Phere*np.exp(-2*1j*phieven)
        betavals.append(np.sum(integrand)/2*np.pi)
    argb2 = -180/np.pi*np.angle(np.array(betavals)) if deg else -np.angle(np.array(betavals))
    return Reven, argb2

def I_interp(im, minR, maxR, minphi, maxphi, narr):
    npix = im.xdim
    fov_muas = im.fovx()/eh.RADPERUAS
    iarr = np.reshape(im.imvec,(npix,npix))
    qarr = np.reshape(im.qvec,(npix,npix))
    uarr = np.reshape(im.uvec,(npix,npix))
    pxi = (np.arange(npix))/npix-0.5
    pxi += (pxi[1]-pxi[0])/2
    pxj = np.copy(pxi)
    mui = pxi*fov_muas
    muj = pxj*fov_muas
    MUI,MUJ = np.meshgrid(mui,muj)
    Reven = np.linspace(minR,maxR,narr)
    phieven = np.arange(minphi,maxphi,2*np.pi/narr)
    
    interp_I = RegularGridInterpolator((mui, muj), iarr)
    interp_Q = RegularGridInterpolator((mui, muj), qarr)
    interp_U = RegularGridInterpolator((mui, muj), uarr)
    
    iintvals = []
    unitcoords = np.array([[-np.cos(phi),-np.sin(phi)] for phi in phieven])
    for R in Reven:
        coords = R*unitcoords
        integrand = interp_I(coords)
        iintvals.append(R*np.sum(integrand)/2*np.pi) #proportional to integrated I in a small annulus
    return Reven, np.array(iintvals)
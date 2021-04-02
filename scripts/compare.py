#!/usr/bin/env python3
"""
  Compare two images and quantify their differences
"""

import numpy as np
import h5py
import sys

# Suppress runtime math warnings -- images have some zeros and that's okay
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")

def blur(im, fwhm=20):
    # 20uas FWHM Gaussian blur. Assume 1px/uas
    return gaussian_filter(im, sigma=(fwhm / (2 * np.sqrt(2 * np.log(2)))))

# Statistics functions
def mse(image1, image2):
    """MSE of image1 vs image2.  As in GRRT paper definition of said, divides by image1 sum"""
    # Avoid dividing by 0
    norm = np.sum(np.abs(image1)**2)
    if norm == 0: norm = 1

    return np.sum(np.abs(image1 - image2)**2) / norm

def ssim(imageI, imageK):
    """SSIM as defined in Gold et al eq 10-ish"""
    N = imageI.shape[0] * imageI.shape[1]
    mu_I = np.sum(imageI) / N
    mu_K = np.sum(imageK) / N
    sig2_I = np.sum((imageI - mu_I)**2) / (N-1)
    sig2_K = np.sum((imageK - mu_K)**2) / (N-1)
    sig_IK = np.sum((imageI - mu_I)*(imageK - mu_K)) / (N-1)
    return (2*mu_I*mu_K) / (mu_I**2 + mu_K**2) * (2*sig_IK) / (sig2_I + sig2_K)

def dssim(imageI, imageK):
    tssim = ssim(imageI, imageK)
    if np.isnan(tssim):
        return 0.0
    else:
        return 1/np.abs(tssim) - 1

if __name__ == "__main__":
    fname1 = sys.argv[1]
    fname2 = sys.argv[2]
    print("Comparing {} vs {}".format(fname1, fname2))

    # load
    with h5py.File(fname1,'r') as hfp:
      unpol1 = np.copy(hfp['unpol']).transpose((1,0))
      imagep = np.copy(hfp['pol']).transpose((1,0,2))
      I1 = imagep[:,:,0]
      Q1 = imagep[:,:,1]
      U1 = imagep[:,:,2]
      V1 = imagep[:,:,3]
    with h5py.File(fname2,'r') as hfp:
      unpol2 = np.copy(hfp['unpol']).transpose((1,0))
      imagep = np.copy(hfp['pol']).transpose((1,0,2))
      I2 = imagep[:,:,0]
      Q2 = imagep[:,:,1]
      U2 = imagep[:,:,2]
      V2 = imagep[:,:,3]

    # command line output
    print("Images\t\tMSE\t\tDSSIM\t\tabs flux diff\trel flux diff")
    mseu = mse(unpol1, unpol2); dssimu = dssim(unpol1, unpol2); diffu = unpol2.sum() - unpol1.sum(); rdiffu = (unpol2.sum() - unpol1.sum())/unpol1.sum()
    print("Unpolarized:\t{:.7g}\t{:.7g}\t{:.7g}\t{:.7g}".format(mseu, dssimu, diffu, rdiffu))
    mseI = mse(I1, I2); dssimI = dssim(I1, I2); diffI = I2.sum() - I1.sum(); rdiffI = (I2.sum() - I1.sum())/I1.sum()
    print("Stokes I:\t{:.7g}\t{:.7g}\t{:.7g}\t{:.7g}".format(mseI, dssimI, diffI, rdiffI))
    mseQ = mse(Q1, Q2); dssimQ = dssim(Q1, Q2); diffQ = Q2.sum() - Q1.sum(); rdiffQ = (Q2.sum() - Q1.sum())/Q1.sum()
    print("Stokes Q:\t{:.7g}\t{:.7g}\t{:.7g}\t{:.7g}".format(mseQ, dssimQ, diffQ, rdiffQ))
    mseU = mse(U1, U2); dssimU = dssim(U1, U2); diffU = U2.sum() - U1.sum(); rdiffU = (U2.sum() - U1.sum())/U1.sum()
    print("Stokes U:\t{:.7g}\t{:.7g}\t{:.7g}\t{:.7g}".format(mseU, dssimU, diffU, rdiffU))
    mseV = mse(V1, V2); dssimV = dssim(V1, V2); diffV = V2.sum() - V1.sum(); rdiffV = (V2.sum() - V1.sum())/V1.sum()
    print("Stokes V:\t{:.7g}\t{:.7g}\t{:.7g}\t{:.7g}".format(mseV, dssimV, diffV, rdiffV))

    lp1 = 100.*np.sqrt(Q1.sum()**2+U1.sum()**2)/I1.sum()
    lp2 = 100.*np.sqrt(Q2.sum()**2+U2.sum()**2)/I2.sum()
    print("Diff LP [%]: {:g}".format(lp2 - lp1))
    cp1 = 100.*V1.sum()/I1.sum()
    cp2 = 100.*V2.sum()/I2.sum()
    print("Diff CP [%]: {:g}".format(cp2 - cp1))
    evpatot1 = 180./3.14159*0.5*np.arctan2(U1.sum(),Q1.sum())
    evpatot2 = 180./3.14159*0.5*np.arctan2(U2.sum(),Q2.sum())
    print("Diff EVPA [deg]: {:g}".format(evpatot2 - evpatot1))

    # Return code for automated testing.  Adjust stringency to taste
    if mseI > 0.005 or mseQ > 0.01 or mseU > 0.01 or mseV > 0.03:
        exit(1)
    else:
        exit(0)

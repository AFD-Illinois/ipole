"""

  Makes images from hdf5 output of ipole.
  2019.07.10 gnw

$ python ipole_plot.py path/to/images/*h5

$ ffmpeg -framerate 8 -i dump%*.png -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4
"""

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys


## configuration / plot parameters

FOV_UNITS = "muas"  # can be set to "muas" or "M" (case sensitive)

## EVPA_CONV not yet implemented! this will only work for observer convention!
EVPA_CONV = "EofN"  # can be set fo "EofN" or "NofW" 



## no need to touch anything below this line

def colorbar(mappable):
  """ the way matplotlib colorbar should have been implemented """
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  ax = mappable.axes
  fig = ax.figure
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  return fig.colorbar(mappable, cax=cax)

if __name__ == "__main__":

  for fname in sys.argv[1:]:

    if fname[-3:] != ".h5": continue
    print("plotting {0:s}".format(fname))

    # load
    hfp = h5py.File(fname,'r')    
    dx = hfp['header']['camera']['dx'][()]
    dsource = hfp['header']['dsource'][()]
    lunit = hfp['header']['units']['L_unit'][()]
    fov_muas = dx / dsource * lunit * 2.06265e11
    scale = hfp['header']['scale'][()]
    evpa_0 = 'W'
    if 'evpa_0' in hfp['header']:
      evpa_0 = hfp['header']['evpa_0'][()]
    unpol = np.copy(hfp['unpol']).transpose((1,0))
    imagep = np.copy(hfp['pol']).transpose((1,0,2))
    I = imagep[:,:,0]
    Q = imagep[:,:,1]
    U = imagep[:,:,2]
    V = imagep[:,:,3]
    hfp.close()

    # set extent (assumption of square image)
    if FOV_UNITS == "muas":
      extent = [ -fov_muas/2, fov_muas/2, -fov_muas/2, fov_muas/2 ]
    elif FOV_UNITS == "M":
      extent = [ -dx/2, dx/2, -dx/2, dx/2 ]
    else:
      print("! unrecognized units for FOV {0:s}. quitting.".format(FOV_UNITS))

    # create plots
    plt.close('all')
    plt.figure(figsize=(8,8))
    ax1 = plt.subplot(2,2,1)
    ax2 = plt.subplot(2,2,2)
    ax3 = plt.subplot(2,2,3)
    ax4 = plt.subplot(2,2,4)

    # get mask for total intensity based on negative values
    Imaskval = np.abs(I.min()) * 100.
    Imaskval = np.nanmax(I) / np.power(I.shape[0],5.)

    # total intensity
    vmax = 1.e-4
    vmax = I.max() / np.sqrt(1.5)
    im1 = ax1.imshow(I, cmap='afmhot', vmin=0., vmax=vmax, origin='lower', extent=extent)
    colorbar(im1)

    # linear polarization fraction
    lpfrac = 100.*np.sqrt(Q*Q+U*U)/I
    lpfrac[np.abs(I)<Imaskval] = np.nan
    ax2.set_facecolor('black')
    im2 = ax2.imshow(lpfrac, cmap='jet', vmin=0., vmax=100., origin='lower', extent=extent)
    colorbar(im2)

    # circular polarization fraction
    cpfrac = 100.*V/I
    cpfrac[np.abs(I)<Imaskval] = np.nan
    vext = max(np.abs(np.nanmin(cpfrac)),np.abs(np.nanmax(cpfrac)))
    vext = max(vext, 1.)
    if np.isnan(vext): vext = 10.
    ax4.set_facecolor('black')
    im4 = ax4.imshow(cpfrac, cmap='seismic', vmin=-vext, vmax=vext, origin='lower', extent=extent)
    colorbar(im4)

    # evpa
    evpa = (180./3.14159)*0.5*np.arctan2(U,Q)
    if evpa_0 == "W":
      evpa += 90.
      evpa[evpa > 90.] -= 180.
    if EVPA_CONV == "NofW":
      evpa += 90.
      evpa[evpa > 90.] -= 180.
    evpa2 = np.copy(evpa)
    evpa2[np.abs(I)<Imaskval] = np.nan
    ax3.set_facecolor('black')
    im3 = ax3.imshow(evpa2, cmap='hsv', vmin=-90., vmax=90., origin='lower', extent=extent)
    colorbar(im3)

    # quiver on intensity
    npix = I.shape[0]
    xs = np.linspace(-fov_muas/2,fov_muas/2,npix)
    Xs,Ys = np.meshgrid(xs,xs)
    lpscal = np.max(np.sqrt(Q*Q+U*U))
    vxp = np.sqrt(Q*Q+U*U)*np.sin(evpa*3.14159/180.)/lpscal
    vyp = -np.sqrt(Q*Q+U*U)*np.cos(evpa*3.14159/180.)/lpscal
    skip = int(npix/32) 
    ax1.quiver(Xs[::skip,::skip],Ys[::skip,::skip],vxp[::skip,::skip],vyp[::skip,::skip], 
      headwidth=1, headlength=1, 
      width=0.005,
      color='#00ff00', 
      units='width', 
      scale=4,
      pivot='mid')

    # command line output
    print("Flux [Jy]:    {0:g} {1:g}".format(I.sum()*scale, unpol.sum()*scale))
    print("I,Q,U,V [Jy]: {0:g} {1:g} {2:g} {3:g}".format(I.sum()*scale,Q.sum()*scale,U.sum()*scale,V.sum()*scale))
    print("LP [%]:       {0:g}".format(100.*np.sqrt(Q.sum()**2+U.sum()**2)/I.sum()))
    print("CP [%]:       {0:g}".format(100.*V.sum()/I.sum()))
    evpatot = 180./3.14159*0.5*np.arctan2(U.sum(),Q.sum())
    if evpa_0 == "W":
      evpatot += 90. 
      if evpatot > 90.:
        evpatot -= 180
    if EVPA_CONV == "NofW":
      evpatot += 90.
      if evpatot > 90.:
        evpatot -= 180
    print("EVPA [deg]:   {0:g}".format(evpatot))

    # formatting and text
    ax1.set_title("Stokes I [cgs]")
    ax2.set_title("LP [%]")
    ax3.set_title("EVPA [deg]")
    ax4.set_title("CP [%]")
    ax1.set_aspect('equal')
    ax2.set_aspect('equal')
    ax3.set_aspect('equal')
    ax4.set_aspect('equal')
    ax1.set_ylabel(FOV_UNITS)
    ax3.set_ylabel(FOV_UNITS)
    ax3.set_xlabel(FOV_UNITS)
    ax4.set_xlabel(FOV_UNITS)

    # saving
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(fname.replace(".h5",".png"))


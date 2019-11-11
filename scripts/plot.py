"""

  Makes images for hdf5 output of ipole.
  last updated: 2018.11.16 gnw

$ python plot.py [-log] [-pol] path/to/images/*h5

$ ffmpeg -framerate 8 -i dump%*.png -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4

"""

print("This is plot.py.")
print("plot.py is decprecated as of 11 Nov 2019.")
print("please use plot_pol.py instead")

import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rc('font', family='serif')
plt.rcParams['mathtext.fontset'] = 'cm'

# if this is not None, then set this as the maximum (linear) value
# for the colorbar. otherwise, set maximum color based on the real
# data.
MAXCOLOR = None

TIME_IN_DAYS = True         # if False, titles plot with time in M
timeoffsetindays = 0.       # if using days, subtract this off
interpolationtype = 'none'  # 'none' or 'bilinear'

COLORMAP = 'afmhot'
DOLOG = False
USE_POLARIZED = False

def colorbar(mappable):
  """ aligns color bar with right of image """
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  ax = mappable.axes
  fig = ax.figure
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  return fig.colorbar(mappable, cax=cax)

if __name__ == "__main__":

  if "-pol" in sys.argv: USE_POLARIZED = True
  if "-log" in sys.argv: DOLOG = True

  fnames = sys.argv[1:]

  for fname in fnames:

    if fname[-3:] != ".h5" and fname[-4:] != ".dat":
      continue

    print("processing {0:s}".format(fname))

    # load data
    hfp = h5py.File(fname,'r') 

    nx = hfp['header']['camera']['nx'][()]
    ny = hfp['header']['camera']['ny'][()]
    dx = hfp['header']['camera']['dx'][()]
    dy = hfp['header']['camera']['dy'][()]
    dsource = hfp['header']['dsource'][()]
    scale = hfp['header']['scale'][()]
    Lunit = hfp['header']['units']['L_unit'][()]
    Munit = hfp['header']['units']['M_unit'][()]
    freqcgs = hfp['header']['freqcgs'][()]

    image = hfp['unpol'][:,:] * scale
    try: tau = hfp['tau'][:,:]
    except: pass
    imageS0 = hfp['pol'][:,:,0] * scale
    imageS1 = hfp['pol'][:,:,1] * scale
    imageS2 = hfp['pol'][:,:,2] * scale
    imageS3 = hfp['pol'][:,:,3] * scale
    imageS4 = hfp['pol'][:,:,4] * scale

    fname = fname.replace(".h5",".dat")

    try:
      timeinm = hfp['header']['t'][()]
    except:
      timeinm = -1
    title = "t = {0:g} M".format(timeinm)

    hfp.close()

    # process/make image and save
    if USE_POLARIZED:
      image = imageS0
    if DOLOG:
      image = np.log10(image)

    plt.close('all')
    fig, (ax) = plt.subplots(1,1)

    fovx_uas = dx * Lunit / dsource * 2.06265e11
    fovy_uas = dy * Lunit / dsource * 2.06265e11
    extent = [ -fovx_uas/2., fovx_uas/2., -fovy_uas/2., fovy_uas/2. ]

    if DOLOG:
      if MAXCOLOR is None: MAXCOLOR = image.max() / 2.
      im = ax.imshow(image.T,interpolation='none',cmap=COLORMAP,origin='lower',vmin=-2,vmax=np.log10(MAXCOLOR))
      plt.title(title)
      colorbar(im)
      plt.tight_layout()
      if USE_POLARIZED:
        plt.savefig(fname.replace('.dat','_log_pol.png'))
      else:
        plt.savefig(fname.replace('.dat','_log.png'))

    else:
      if MAXCOLOR is None: MAXCOLOR = image.max()

      # note that this rotation makes it so the jet comes out to the right
      im = ax.imshow(image.T,interpolation=interpolationtype,cmap=COLORMAP,origin='lower',vmax=MAXCOLOR,extent=extent)

      if timeinm >= 0:
        if TIME_IN_DAYS: # plots time in days
          timeinhours = timeinm * Lunit / 2.998e10 / 3600.
          timeindays = timeinhours / 24 - timeoffsetindays
          ttext = "$\mathrm{+\," + "{0:.01f}".format(np.abs(timeindays)) + "\ days}$"
          plt.suptitle(ttext,fontsize=17)
        else: # plots time in M
          ttext = "$\mathrm{t = " + str(int(timeinm)) + "\ M}$"
          plt.suptitle(ttext,fontsize=17)

      plt.xlabel(r"$\mathrm{x}\;\;[\mu{}\mathrm{as}]$",fontsize=15)
      plt.ylabel(r"$\mathrm{y}\;\;[\mu\mathrm{as}]$",fontsize=15)
      plt.tight_layout(rect=[0,0.03,1,0.95])

      if USE_POLARIZED:
        plt.savefig(fname.replace('.dat','_pol.png'))
      else:
        plt.savefig(fname.replace('.dat','.png'),dpi=200)


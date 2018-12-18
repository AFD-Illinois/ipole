"""
  
  Utility that can be used to translate from native ipole output format to 
  FITS image file format for EHT theory/model group. Generally descended
  from GRRT2fits script provided by C.M. Fromm.

  last modified: 2018.11.17 gnw

  $ python translate-fits.py path/to/ipole/files/*.h5

"""

# general configuration parameters
user = "C.F. Gammie"
data_to_save = "unpolarized"  # "polarized" or "unpolarized"
source = "M87"                # "M87" or "SgrA*"
date = "2018-11-16"
code = "iharm3d"


# no need to touch anything below this line
import sys
import numpy as np
import h5py
from astropy.time import Time as aTime
import astropy.io.fits as fits

def get_small_float_string(num):
  """ returns a string with no trailing zeros """
  return "{0:f}".format(num).rstrip('0').rstrip('.')

def write_fits_from_ipole(fname, ofname):
  """ translates ipole file at 'fname' to fits file at 'ofname' """

  hfp = h5py.File(fname,'r')

  # set various configuration attributes
  ra = None
  dec = None
  if source == "M87" and True:
    ra = 187.70593075
    dec = 12.391123306
  else:
    from astropy.coordinates import SkyCoord
    loc = SkyCoord.from_name(source)
    ra = loc.ra.deg
    dec = loc.dec.deg
  DEGREE = 3.141592653589/180.0
  HOUR = 15.0*DEGREE
  RADPERAS = DEGREE/3600.0
  RADPERUAS = RADPERAS*1.e-6
  modjuldate=(aTime(date)).mjd

  # read header data from file
  t = hfp['header']['t'][()]
  freq = hfp['header']['freqcgs'][()]
  scale = hfp['header']['scale'][()]
  nx = hfp['header']['camera']['nx'][()]
  ny = hfp['header']['camera']['ny'][()]
  DX = hfp['header']['camera']['dx'][()]
  DY = hfp['header']['camera']['dx'][()]
  rcam = hfp['header']['camera']['rcam'][()]
  thetacam = hfp['header']['camera']['thetacam'][()]
  phicam = hfp['header']['camera']['phicam'][()]
  Lunit = hfp['header']['units']['L_unit'][()]
  Dsource = hfp['header']['dsource'][()]
  fovx_uas = DX * Lunit / Dsource * 2.06265e11
  fovy_uas = DY * Lunit / Dsource * 2.06265e11
  dxorg = -1. * fovx_uas / nx
  dyorg = 1. * fovy_uas / ny

  fheader = hfp['fluid_header']
  metric = fheader['metric'][()].lower()
  a = fheader['geom'][metric]['a'][()]
  inc = get_small_float_string(thetacam)
  time = "{0:d}M".format(int(t))
  spin = get_small_float_string(a)
  resolution = "{0:d}x{1:d}x{2:d}".format(fheader['n1'][()],fheader['n2'][()],fheader['n3'][()])

  if hfp['header']['electrons']['type'][()] == 0:
    hstring = "tp/te={0:s}".format(get_small_float_string(hfp['header']['electrons']['tp_over_te'][()]))
  elif hfp['header']['electrons']['type'][()] == 1:
    hstring = "(kawazura Te)"
  elif hfp['header']['electrons']['type'][()] == 2:
    hstring = "Rlow={0:s},Rhigh={1:s}".format(\
      get_small_float_string(hfp['header']['electrons']['rlow'][()]),\
      get_small_float_string(hfp['header']['electrons']['rhigh'][()]))
  else:
    hstring = "(unlisted Te)"

  history = """GRMHD: {0:s}, Kerr, a={1:s}, 3D, {2:s}""".format(code, spin, resolution)
  history += """ GRRT: ipole, t={0:s}, {1:s}, i={2:s}""".format(time, hstring, inc)

  # read flux data from file, rotate according to spin to aim the jet right
  unpol_Jy = hfp['unpol'][:,:].T * scale
  I_Jy = hfp['pol'][:,:,0].T * scale

  if "un" in data_to_save:
    data = unpol_Jy
  else:
    data = I_Jy

  if a >= 0:
    data = np.rot90(data,k=3)
  else:
    data = np.rot90(data,k=1)

  hfp.close()

  # all code below this line is associated with actually writing the
  # fits file

  # fill out the header
  header = fits.Header()
  header['AUTHOR'] = user
  header['OBJECT'] = source
  header['CTYPE1'] = 'RA---SIN'
  header['CTYPE2'] = 'DEC--SIN'
  header['CDELT1'] = dxorg*RADPERUAS/DEGREE
  header['CDELT2'] = dyorg*RADPERUAS/DEGREE
  header['OBSRA'] = ra
  header['OBSDEC'] = dec
  header['FREQ'] = freq
  header['MJD'] = float(modjuldate)
  header['TELESCOP'] = 'VLBI'
  header['BUNIT'] = 'JY/PIXEL'
  header['STOKES'] = 'I'
  header['HISTORY'] = history

  # create file with data
  hdu = fits.PrimaryHDU(data, header=header)
  hdulist = [hdu]
  hdulist = fits.HDUList(hdulist)

  # save the file
  hdulist.writeto(ofname, overwrite=True)


if __name__ == "__main__":

  fnames = [ x for x in sys.argv[1:] if x[-3:] == ".h5" ]

  for fname in fnames:

    print("translating {0:s}".format(fname))

    ofname = fname.replace(".h5",".fits")  # TODO final format of FITS filename should be 'FITS/%s_%i_%iGHz.fits' %(outname[:-4],dim,(freq/1e9))
    write_fits_from_ipole(fname,ofname)

    

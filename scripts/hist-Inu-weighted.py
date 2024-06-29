"""

  You probably want to use this file ...

  2019.01.16

"""

import os
import sys
import h5py
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pyharm

# set parameters here
rlow = 1.
rhigh = 80.
thetae_mode = "MIXED" # "FIXED" "MIXED" "NATIVE"  (only MIXED implemented right now)
sinth = 1./2


# no need to touch anything below this line

QE = 4.80320680e-10  # cgs
CL = 2.99792458e10   # cgs
ME = 9.1093826e-28   # cgs
MP = 1.67262171e-24  # cgs

def plthist(title,fname,data,wgts,dolog=True,nbins=100,vmin=None,vmax=None,ymin=None,ymax=None):
  plt.close('all')
  if vmin is None: vmin = data.min()
  if vmax is None: vmax = data.max()
  plt.hist(data,weights=wgts,log=dolog,bins=np.logspace(np.log10(vmin),np.log10(vmax),nbins))
  plt.xscale('log')
  if ymin is not None: plt.ylim(bottom=ymin)
  if ymax is not None: plt.ylim(top=ymax)
  plt.title(title)
  plt.savefig(fname)

def colorbar(mappable):
  """ the way matplotlib colorbar should have been implemented """
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  ax = mappable.axes
  fig = ax.figure
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  return fig.colorbar(mappable, cax=cax)

def get_ucon(U, gcon, gcov):
  N1 = U.shape[0]
  N2 = U.shape[1]
  N3 = U.shape[2]
  alpha = 1. / np.sqrt(-gcon[:,:,:,0,0])
  gamma = np.sqrt(1. + np.einsum('abci,abci->abc',np.einsum('abcij,abci->abcj',gcov[:,:,:,1:,1:],U),U))
  ucon = np.zeros((N1,N2,N3,4))
  ucon[:,:,:,1:] = U - gamma[:,:,:,None] * alpha[:,:,:,None] * gcon[:,:,:,0,1:]
  ucon[:,:,:,0] = gamma / alpha
  return ucon

if __name__ == "__main__":

  # set these here
  dump_fname = "tests/test-resources/sample_dump_SANE_a+0.94_MKS_0900.h5"
  zone_fname = "histo.h5"
  imag_fname = "image.h5"

  # load grid values (for computation of four-vector quantities)
  dump = pyharm.load_dump(dump_fname)

  # load from dump file
  print("loading dump file now ...")
  hfp = h5py.File(dump_fname,'r')
  N1 = hfp["header"]["n1"][()]
  N2 = hfp["header"]["n2"][()]
  N3 = hfp["header"]["n3"][()]
  startx1 = hfp["header"]["geom"]["startx1"][()]
  dx1 = hfp["header"]["geom"]["dx1"][()]
  gam = hfp["header"]["gam"][()]
  game = 4./3
  gamp = 5./3
  hfp.close()
  (RHO, UU, U, B, ucon, ucov, bcon, bcov) = dump['rho'],  dump['u'],  dump['uvec'],  dump['Bvec'],  dump['ucon'],  dump['ucov'],  dump['bcon'],  dump['bcov']
  print("done!")

  # get units
  with h5py.File(imag_fname,'r') as imag:
    Lunit = imag["header"]["units"]["L_unit"][()]
    Munit = imag["header"]["units"]["M_unit"][()]
    Tunit = imag["header"]["units"]["T_unit"][()]
    scale = imag["header"]["scale"][()]

  RHOunit = Munit/np.power(Lunit,3.)
  Bunit = CL * np.sqrt(4.*np.pi*RHOunit)

  # calculate physical quantities
  BSQ = np.einsum('iabc,iabc->abc',bcon,bcov)
  BETA = UU*(gam-1.)/0.5/BSQ
  BETASQ = np.power(BETA,2.)
  TRAT = rhigh*BETASQ/(1.+BETASQ) + rlow*BETASQ/(1.+BETASQ)
  THETAE = None
  if thetae_mode == "MIXED":
    THETAE = UU/RHO * MP/ME * (game-1.)*(gamp-1.)/((gamp-1.)+(game-1.)*TRAT)
  NUS = QE * np.sqrt(BSQ)*Bunit * THETAE*THETAE * sinth / (9.*np.pi*ME*CL)
  SYNCH_X = 230.e+9 / NUS

  # load from zone file
  with h5py.File(zone_fname,'r') as zones:
    zdata = zones["data"][()] * np.power(230.e9, 3.) * scale

  # generate geometry
  X1i,X2i,X3i = np.meshgrid(range(N1),range(N2),range(N3))
  X1i = X1i.transpose((1,0,2))
  X2i = X2i.transpose((1,0,2))
  X3i = X3i.transpose((1,0,2))
  r = np.exp(np.linspace(startx1,startx1+dx1*N1,N1))
  h = np.linspace(0,np.pi,N2) # TODO: make this right
  p = np.linspace(0,2.*np.pi,N3)
  R,H,P = np.meshgrid(r,h,p)
  R = R.transpose((1,0,2))
  H = H.transpose((1,0,2))
  P = P.transpose((1,0,2))

  os.makedirs("imgs", exist_ok=True)
  # r
  plthist("R [M]","imgs/r-hist.png",R.flatten(),zdata.flatten(),ymin=1.e-4,ymax=1,vmin=2,vmax=100)

  # B
  plthist("B [cgs]","imgs/b-hist.png",np.sqrt(BSQ.flatten())*Bunit,zdata.flatten(),vmin=1.e-2,vmax=1.e+2,ymin=1.e-10)

  # rho
  plthist("RHO [cgs]","imgs/rho-hist.png",RHO.flatten()*RHOunit,zdata.flatten())

  # Thetae
  plthist("THETAE","imgs/thetae-hist.png",THETAE.flatten(),zdata.flatten(),vmin=1.e-3,vmax=1.e4)
  
  # x = nu / nu_s
  # nus = e B Thetae^2 Sin[th] / (9 Pi me c)
  plthist("nu/nus","imgs/x-hist.png",SYNCH_X.flatten(),zdata.flatten(),vmin=1.e-5,vmax=1.e+10)

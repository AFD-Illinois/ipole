from subprocess import call
import numpy as np
import pickle
from glob import glob
from scipy.ndimage.interpolation import rotate

CL = 2.99792458e10

theta = 30.
nu = 230.e9
phi = 288.
DTd = 5.
folder = '/data/bh-fs4/bryan10/m87_2d/M3e9/a09/dumps/'

files = np.sort(glob(folder+'dump_fluid*'))

times = np.zeros(len(files))
nuFnu = np.zeros(len(files))
nuFnu_unpol = np.zeros(len(files))
stokes = []

for n, fnam in enumerate(files):
  print n
  print fnam
  call(['./ipole', str(theta), str(nu), fnam, '1', '1', '1'])
  with open('ipole.dat', 'r') as f:
    line = f.readline().split(' ')
    N1 = int(line[0])
    N2 = int(line[1])
    DX = float(line[2])
    DY = float(line[3])
    scale = float(line[4])
    L_unit = float(line[5])
    M_unit = float(line[6])
    
  i0, j0, x, y, Ia, Is, Qs, Us, Vs = np.loadtxt('ipole.dat', unpack=True,
      skiprows=1)
  N = int(np.sqrt(len(i0)))
  N1 = N2 = N
 
  times[n] = n*DTd
  nuFnu[n] = np.sum(Is)*scale
  nuFnu_unpol[n] = np.sum(Ia)*scale
  print nuFnu[n]
  
  x = np.reshape(x, (N1, N2)) - 20
  y = np.reshape(y, (N2, N2)) - 20
  Is = rotate(np.reshape(Is, (N1, N2)), phi, reshape=False)
  Qs = rotate(np.reshape(Qs, (N1, N2)), phi, reshape=False)
  Us = rotate(np.reshape(Us, (N1, N2)), phi, reshape=False)
  Vs = rotate(np.reshape(Vs, (N1, N2)), phi, reshape=False)
  stokes.append([x, y, Is, Qs, Us, Vs])

out = {}

out['nuFnu [Jy]'] = nuFnu
out['nuFnu_unpol [Jy]'] = nuFnu_unpol
#out['Stokes'] = np.array(stokes)
out['N'] = N
out['nu'] = nu
out['theta'] = theta
out['phi'] = phi
out['DTd'] = DTd
out['t'] = times
out['folder'] = folder

pickle.dump(out, open('m87.p', 'wb'))


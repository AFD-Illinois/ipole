from subprocess import call
import numpy as np
import pickle
from glob import glob
from scipy.ndimage.interpolation import rotate

CL = 2.99792458e10

theta = 20.
nu = 230.e9
phi = 288.
DTd = 5.
mass = 'M3e9'
spin = 'a09'
folder = '/data/bh-fs4/bryan10/m87_2d/' + mass + '/' + spin + '/dumps/'

files = np.sort(glob(folder+'dump_fluid*'))
files = files[120:]

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
    NX = int(line[0])
    NY = int(line[1])
    DX = float(line[2])
    DY = float(line[3])
    scale = float(line[4])
    L_unit = float(line[5])
    M_unit = float(line[6])
    
  i0, j0, x, y, Ia, Is, Qs, Us, Vs, tauF = np.loadtxt('ipole.dat', unpack=True,
      skiprows=1)
  #N = int(np.sqrt(len(i0)))
  #NX = NY = N
 
  times[n] = n*DTd
  nuFnu[n] = np.sum(Is)*scale
  nuFnu_unpol[n] = np.sum(Ia)*scale
  print nuFnu[n]
  
  x = np.reshape(x, (NX, NY)) - DX/2
  y = np.reshape(y, (NX, NY)) - DY/2
  Is = rotate(np.reshape(Is, (NX, NY)), phi, reshape=False)
  Qs = rotate(np.reshape(Qs, (NX, NY)), phi, reshape=False)
  Us = rotate(np.reshape(Us, (NX, NY)), phi, reshape=False)
  Vs = rotate(np.reshape(Vs, (NX, NY)), phi, reshape=False)
  stokes.append([x, y, Is, Qs, Us, Vs])

out = {}

out['nuFnu [Jy]'] = nuFnu
out['nuFnu_unpol [Jy]'] = nuFnu_unpol
#out['Stokes'] = np.array(stokes)
out['N'] = NX
out['nu'] = nu
out['theta'] = theta
out['phi'] = phi
out['DTd'] = DTd
out['t'] = times
out['folder'] = folder

pickle.dump(out, open('m87_' + mass + '_' + spin + '.p', 'wb'))


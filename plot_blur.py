import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
  print('ERROR format is')
  print('  python plot.py [ipole.dat]')
  sys.exit()

fnam = sys.argv[1]
print(fnam)

data = h5py.File(fnam, 'r')
print(list(data.keys()))
print(data['nuLnu'][()])

DX = np.copy(data['/header/camera/dx'])
NX = np.copy(data['/header/camera/nx'])
D = np.copy(data['/header/dsource'])
L_unit = np.copy(data['/header/units/L_unit'])
#print(DX, NX)

uas_per_pix = DX*L_unit/NX/D/4.8481368110954E-12
M_to_uas = L_unit/D/4.8481368110954E-12
print(DX,  D)
#print(uas_per_pix)
print(M_to_uas)

#sys.exit()

unpol = np.copy(data['unpol'])
pol = np.copy(data['pol'])
X = np.copy(data['X'])
Y = np.copy(data['Y'])
data.close()

plt.figure()
ax = plt.gca()

from scipy.ndimage.interpolation import rotate
Z = rotate(unpol, 288, reshape=False)

ax = plt.subplot(1,2,1)
ax.pcolormesh(X, Y, Z, cmap='afmhot')
ax.set_aspect('equal')


ax = plt.subplot(1,2,2)
from scipy.ndimage.filters import gaussian_filter as gf
Z = gf(Z, 15/uas_per_pix/np.sqrt(8*np.log(2.)), truncate=3.0)
ax.pcolormesh(X*M_to_uas, Y*M_to_uas, Z, cmap='afmhot')

#ax.pcolormesh(X, Y, pol[:,:,0], cmap='afmhot')
ax.set_aspect('equal')
plt.show()



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
unpol = np.copy(data['unpol'])
pol = np.copy(data['pol'])
X = np.copy(data['X'])
Y = np.copy(data['Y'])
data.close()

plt.figure()
ax = plt.gca()
ax.pcolormesh(X, Y, unpol, cmap='afmhot')
#ax.pcolormesh(X, Y, pol[:,:,0], cmap='afmhot')
ax.set_aspect('equal')
plt.show()



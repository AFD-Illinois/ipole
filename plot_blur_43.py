import h5py
import numpy as np
import sys
import matplotlib
font = {'family' : 'serif',
        'size'   : 18}
matplotlib.rc('font', **font)
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

#fig = plt.figure(figsize=(14,10))
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14,4))
#ax = plt.gca()

from scipy.ndimage.interpolation import rotate
Z = rotate(unpol, 288-180, reshape=False)
Z = np.clip(Z, 1.e-10*Z.max(), a_max = None)

ax = axes[0]
#ax = plt.subplot(1,2,1)
ax.pcolormesh(X*M_to_uas, Y*M_to_uas, np.log10(Z/Z.max()), cmap='jet', vmin=-4, vmax=0)
ax.set_aspect('equal')
ax.set_xlim([-1000,1000])
ax.set_ylim([-1000,1000])
ax.set_xlabel('X [uas]')
ax.set_ylabel('Y [uas]')


ax = axes[1]
#ax = plt.subplot(1,2,2)
from scipy.ndimage.filters import gaussian_filter as gf
sigma = [210/2/uas_per_pix/np.sqrt(8*np.log(2.)), 430/2/uas_per_pix/np.sqrt(8*np.log(2.))]
Z = rotate(Z, 16, reshape=False)
Z = gf(Z, sigma, truncate=3.0)
Z = rotate(Z, -16, reshape=False)
Z = np.log10(Z/Z.max())
cax = ax.pcolormesh(X*M_to_uas, Y*M_to_uas, Z, cmap='jet', vmin=-4, vmax=0)
fig.colorbar(cax, ax=axes.ravel().tolist(), label='log Intensity [arb]')
ax.contour(X*M_to_uas, Y*M_to_uas, Z, colors=['w'], levels=[-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0],
  linestyles=['-'], linewidths=[1])
ax.set_xlim([-1000,1000])
ax.set_ylim([-1000,1000])
ax.set_yticklabels([])
ax.set_xlabel('X [uas]')

#ax.pcolormesh(X, Y, pol[:,:,0], cmap='afmhot')
ax.set_aspect('equal')
plt.show()
#plt.savefig('/home/brryan/Documents/EHT/nov2018_conf_figs/M87_43GHz.png', dpi=300, bbox_inches='tight')



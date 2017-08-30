from subprocess import call
import numpy as np
import pickle
import matplotlib.pyplot as plt
from glob import glob

CL = 2.99792458e10

theta = 20.
nu = 230.e9
N = 128
folder = '/data/bh-fs4/bryan10/m87_2d/M3e9/a05/dumps/'

#fnam = 'dump_fluid00000200'
#EVPA = np.zeros(len(nu))

files = np.sort(glob(folder+'dump_fluid*'))

nuFnu = np.zeros(len(files))

for n, fnam in enumerate(files):
  print n
  print fnam
  call(['./ipole', str(theta), str(nu), fnam, '1', '1', '1'])
  i0, j0, x, y, Ia, Is, Qs, Us, Vs = np.loadtxt('ipole.dat', unpack=True)
  nuFnu[n] = np.sum(Is)

#for n in xrange(len(nu)):
  #call(['./ipole', str(theta), str(nu[n]), fnam, '1', '1', '1'])
  
  #i0, j0, x, y, Ia, Is, Qs, Us, Vs = np.loadtxt('ipole.dat', unpack=True)
  #Q_I = np.sum(Qs*Is)/np.sum(Is)
  #U_I = np.sum(Us*Is)/np.sum(Is)
  #EVPA[n] = 0.5*np.arctan(U_I/Q_I)
  
  #call(['mv', 'ipole.dat', 'ipole_' + '%08d' % n + '.dat'])
  #print nu
  print nuFnu[n]

out = {}
out['nu [Hz]'] = nu
out['lam [cm]'] = CL/nu
out['lam [m'] = out['lam [cm]']*0.01
out['EVPA [rad]'] = EVPA
out['EVPA [deg]'] = EVPA*180./np.pi

print EVPA

pickle.dump(out, open('evpa.p', 'wb'))

ax = plt.subplot(1,1,1)
plt.plot(out['lam [m']**2, out['EVPA [rad]'], marker='s', linestyle='',
  color='k')
ax.set_xlabel('lam^2 (m^2)'); ax.set_ylabel('EVPA (rad)')
ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))

plt.show()


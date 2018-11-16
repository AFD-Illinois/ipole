from subprocess import call
import numpy as np
import pickle
import matplotlib.pyplot as plt

CL = 2.99792458e10

theta = 20.
#nu = np.array([229.6e9, 229.7e9, 229.8e9, 229.9e9, 230.e9, 230.1e9, 230.2e9, 230.3e9, 230.4e9])
nu = np.linspace(35.e9, 700.e9, 20)
nu = np.logspace(np.log10(35.e9), np.log10(700.e9), 128)
fnam = 'dump_fluid00000200'
EVPA = np.zeros(len(nu))

for n in xrange(len(nu)):
  call(['./ipole', str(theta), str(nu[n]), fnam, '1', '1', '1'])
  
  i0, j0, x, y, Ia, Is, Qs, Us, Vs, tauF = np.loadtxt('ipole.dat', unpack=True,
      skiprows=1)
  Q_I = np.sum(Qs*Is)/np.sum(Is)
  U_I = np.sum(Us*Is)/np.sum(Is)
  EVPA[n] = 0.5*np.arctan(U_I/Q_I)
  
  call(['mv', 'ipole.dat', 'ipole_' + '%08d' % n + '.dat'])
  print nu

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


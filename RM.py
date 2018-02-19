from subprocess import call
import numpy as np
import pickle
from glob import glob
from scipy.ndimage.interpolation import rotate
import sys

CL = 2.99792458e10

if len(sys.argv) != 4:
  print 'Format is python RM_counterjet.py [nu] [theta] [dumpfile]'
  exit(-1)

nu    = float(sys.argv[1])
theta = float(sys.argv[2])
dfnam = sys.argv[3]

# Calculate rotation measure
#nus = np.array([0.998*nu, 0.999*nu, nu, 1.001*nu, 1.002*nu])
#nus = np.array([218.4e9, 220.4e9, 230.3e9, 232.3e9])
nus = np.array([228,229,230,231,232])*1.e9
#nus = np.linspace(220.e9, 230.e9, 16)
lams = CL/nus/100. # cm -> m
chis = np.zeros(len(nus))
for n, freq in enumerate(nus):
  print '\n\n\n%i' % n
  print 'freq = %e' % freq
  call(['./ipole', str(theta), str(freq), dfnam, '1', '1', '1', '0'])
  i0, j0, x, y, Ia, Is, Qs, Us, Vs, tauF = np.loadtxt('ipole.dat', unpack=True,
    skiprows=1)
  print '\n\n\n\n'
  N = int(np.sqrt(len(i0)))
  Q = sum(Qs)#sum(Qs*Is)/sum(Is)
  U = sum(Us)#sum(Us*Is)/sum(Is)
  print N
  print Q
  print U
  print 'I = %e' % sum(Is)
  chis[n] = 0.5*np.arctan(U/Q)

RM = np.polyfit(lams**2, chis, 1)[0]
print 'N = %d' % N
print 'RM = %g rad m^-2' % RM
print (lams**2)
print chis

RM_out = {}
RM_out['nu'] = nus
RM_out['lam2'] = lams**2
RM_out['chis'] = chis
RM_out['RM'] = RM
RM_out['dfnam'] = dfnam
RM_out['theta'] = theta
pickle.dump(RM_out, open('RM.p', 'wb'))


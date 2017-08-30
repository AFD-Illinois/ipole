from subprocess import call
import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
  print 'Format: python batch_plot.py [filename]'
  sys.exit()

fname = sys.argv[1]

#fname = 'evpa.p'
data = pickle.load(open(fname, 'rb'))

#['nu [Hz]'] = nu
#['lam [cm]'] = CL/nu
#['lam [m'] = out['lam [cm]']*0.01
#['EVPA [rad]'] = EVPA
#['EVPA [deg]'] = EVPA*180./np.pi

# Update phase
EVPA = np.copy(data['EVPA [rad]'])
phase = 0.
for n in range(1, len(EVPA)):
  if data['EVPA [rad]'][n-1] > np.pi/8. and data['EVPA [rad]'][n] < -np.pi/8.:
    phase += np.pi/2.
  EVPA[n] += phase

ax = plt.subplot(2,1,1)
#plt.plot(data['lam [m']**2, data['EVPA [rad]'], marker='s', linestyle='',
#  color='k')
#ax.set_xlabel('lam^2 (m^2)'); ax.set_ylabel('EVPA (rad)')
#plt.plot(data['nu [Hz]'], data['EVPA [rad]'], marker='s', linestyle='',
#  color='k')
plt.plot(data['nu [Hz]'], EVPA, linestyle='-',
  color='k')
#plt.plot(data['nu [Hz]'], data['EVPA [rad]'], linestyle='-',
#  color='k')
ax.set_xlabel('nu (Hz)'); ax.set_ylabel('EVPA (rad)')
ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))
ax.set_xscale('log')

ax = plt.subplot(2,1,2)
# Get numerical derivative
#EVPA = data['EVPA [rad]']
lam2 = data['lam [m']**2
deriv = np.zeros(len(EVPA))
for n in range(1,len(EVPA) - 1):
  deriv[n] = (EVPA[n+1] - EVPA[n-1])/(lam2[n+1] - lam2[n-1])
plt.plot(data['nu [Hz]'], deriv, color='k')
ax.set_xlabel('nu (Hz)'); ax.set_ylabel('RM (rad m^-2)')
ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))
ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))
ax.set_xscale('log')


plt.show()


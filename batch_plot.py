from subprocess import call
import numpy as np
import pickle
import matplotlib.pyplot as plt

fname = 'evpa.p'
data = pickle.load(open(fname, 'rb'))

#['nu [Hz]'] = nu
#['lam [cm]'] = CL/nu
#['lam [m'] = out['lam [cm]']*0.01
#['EVPA [rad]'] = EVPA
#['EVPA [deg]'] = EVPA*180./np.pi

ax = plt.subplot(1,1,1)
plt.plot(data['lam [m']**2, data['EVPA [rad]'], marker='s', linestyle='',
  color='k')
ax.set_xlabel('lam^2 (m^2)'); ax.set_ylabel('EVPA (rad)')
ax.get_xaxis().get_major_formatter().set_powerlimits((0, 0))
ax.set_xscale('log')

plt.show()


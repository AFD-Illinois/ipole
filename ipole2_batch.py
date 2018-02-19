import sys
from subprocess import call
import numpy as np
import glob

# Loop over time
folder = '/data/bh-fs4/bryan10/m87_2d/M3e9/a09/dumps/'
files = np.sort(glob.glob(folder+'dump_fluid*'))
files = files[120:-1]

#for n, fnam in enumerate(files):
#  call(['./ipole','20','230.e9',fnam,'1','1','1'])
#  call(['python','ipole2.py','ipole.dat','%d' % n])

# Loop over frequency
nus = np.logspace(np.log10(30.e9), np.log10(700.e9), 80)
fnam = files[0]
for n, nu in enumerate(nus):
  call(['./ipole','20',str(nu),fnam,'1','1','1'])
  call(['python','ipole2_nu.py','ipole.dat','%d' % n, str(nu)])


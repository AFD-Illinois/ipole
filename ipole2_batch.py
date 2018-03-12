import sys
from subprocess import call
import numpy as np
  
# Loop over frequency
nus = np.logspace(np.log10(30.e9), np.log10(700.e9), 300)
fnam = '/data/bh-bd2/bryan10/grmhd/premad/dumps/dump_00008000.h5'
for n, nu in enumerate(nus):
  call(['./ipole','80',str(nu),fnam,'5.e16','0','0','0'])
  call(['mv','ipole.dat','ipole_%08d.dat' % n])


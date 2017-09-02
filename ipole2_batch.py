import sys
from subprocess import call
import numpy as np
import glob

folder = '/data/bh-fs4/bryan10/m87_2d/M3e9/a09/dumps/'
files = np.sort(glob.glob(folder+'dump_fluid*'))

print files

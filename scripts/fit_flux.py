#!/usr/bin/env python3

# Use: ./fit_flux.py parfile.par target_flux guess_munit

import sys
import numpy as np

sys.path.append("scripts")
import ipole

parfile = sys.argv[1]
target = float(sys.argv[2])
start = float(sys.argv[3])

run_unpol = lambda munit: ipole.run({'M_unit': munit}, unpol=True, quench=True, parfile=parfile, verbose=1)['Ftot_unpol'] - target

x1 = start
x2 = 10*x1

y1 = run_unpol(x1)
y2 = run_unpol(x2)

while np.abs(y2) > 0.0001:
    xnew = x2 - y2 * (x2 - x1) / (y2 - y1)
    x1 = x2
    x2 = xnew
    y1 = y2
    y2 = run_unpol(x2)
    print("M_unit: {} Flux: {}".format(x2, y2+target))

ipole.run({'M_unit': x2}, parfile=parfile)

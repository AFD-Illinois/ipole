# Run the (unpolarized) GRRT test problems and compare to published values

import os, sys
import numpy as np
sys.path.append('../../scripts/')
from ipole import run

published = [0, 1.6465, 1.4360, 0.4418, 0.2710, 0.0255]

bad = 0
for i in range(1,6):
    results = run({'model': i}, parfile="analytic.par")
    err = np.abs(results['Ftot_unpol'] / published[i] - 1)
    print("Ftot: {} Ftot_pol: {} Exact value: {} Relative Error: {}".format(results['Ftot_unpol'], results['Ftot_pol'], published[i], err))
    if err > 0.02:
        bad = 1

exit(bad)

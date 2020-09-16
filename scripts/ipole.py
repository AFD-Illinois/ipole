"""

  ipole.py (requires python3)

  provides 
   + hooks into running ipole 
   + (TODO) method to read (more of) ipole output to python dictionary

  2018.10.30 gnw

"""

import subprocess

def run(args, exe="./ipole", quench=False, unpol=False, parfile=None, verbose=0):
  """Runs ipole with config as specified by args."""

  cmd = [exe]
  if parfile is not None:
    cmd += ["-par",parfile]

  cmd += ["--{}={}".format(key,args[key]) for key in args]

  if quench: cmd += ["-quench"]
  if unpol: cmd += ["-unpol"]

  if verbose>0: print(" ".join(cmd))
  proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output = [ z for y in [ str(x)[2:-1].split("\\n") for x in proc.communicate() ] for z in y ]

  results = {}
  for line in output:
    if verbose>1: print(line)
    if "Ftot" in line:
      proc = line.replace('(','').replace(')','').split()
      results['Ftot_pol'] = float(proc[3])
      results['Ftot_unpol'] = float(proc[5])

  return results

def run_legacy(thetacam, freqcgs, Mbh, Munit, fname, Rlow=None, Rhigh=None, exe="./ipole", counterjet=0, 
        quench=False, verbose=False, unpol=False):
  """ runs ipole with config as specified by input arguments """
  if Rlow is None and Rhigh is not None: Rlow = 1.
  if Rhigh is None and Rlow is not None: Rhigh = 1.
  if Rlow is None:
    args = [ exe, thetacam, freqcgs, Mbh, Munit, fname, counterjet ]
  else:
    args = [ exe, thetacam, freqcgs, Mbh, Munit, fname, counterjet, Rlow, Rhigh ]
  args = [ str(x) for x in args ]
  if quench: args.append("-quench")
  if unpol: args.append("-unpol")
  if verbose: print(args)
  proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output = [ z for y in [ str(x)[2:-1].split("\\n") for x in proc.communicate() ] for z in y ]
  results = {}
  for line in output:
    if verbose: print(line)
    if "Ftot" in line:
      proc = line.replace('(','').replace(')','').split()
      results['Ftot_pol'] = float(proc[3])
      results['Ftot_unpol'] = float(proc[4])
  return results


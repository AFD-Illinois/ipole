################################################################################
#                                                                              #
#  CONFIGURATION AND COMPILATION ROUTINE                                       #
#                                                                              #
################################################################################

import sys; sys.dont_write_bytecode=True
import machines
import os
import numpy as np
import glob
from subprocess import call
import subprocess

DEBUG = 0
EXECUTABLE = 'ipole'

class color:
  BOLD    = '\033[1m'
  WARNING = '\033[1;31m'
  BLUE    = '\033[94m'
  NORMAL  = '\033[0m'

def build():
  NOPARAM = 1
  DEBUG = 0
  for n in range(len(sys.argv)):
    if sys.argv[n] == '-noparam':
      NOPARAM = 1
    if sys.argv[n] == '-debug':
      DEBUG = 1

  print("")
  print("********************************************************************************")
  print("")
  print("                                  BUILD SCRIPT")
  print("")
  print("********************************************************************************")
  print("")

  host = machines.get_machine()

  C_FLAGS = '-std=c99 ' + host['COMPILER_FLAGS']

  # MATH AND DYNAMIC LINKING
  LIB_FLAGS = '-lm -ldl'

  # GSL
  host['GSL_DIR'] = os.path.join(host['GSL_DIR'], '')
  LIB_FLAGS += ' -lgsl -lgslcblas'
  LIBRARIES = '-L' + host['GSL_DIR'] + 'lib/'
  INCLUDES  = '-I' + host['GSL_DIR'] + 'include/'

  print("  CONFIGURATION\n")
  def print_config(key, var):
    print("    " + color.BOLD + "{:<15}".format(key) + color.NORMAL +
          str(var))

  print_config("MACHINE",    host['NAME'])
  print_config("COMPILER",   host['COMPILER'])
  print_config("GSL_DIR",    host['GSL_DIR'])
  print_config("C_FLAGS",    C_FLAGS)
  print_config("LIB_FLAGS",  LIB_FLAGS)
  print_config("LIBRARIES",  LIBRARIES)
  print_config("INCLUDES",   INCLUDES)

  # GET SINGLE LISTS OF ALL SOURCE, OBJECT, AND HEADER FILES
  SRC_ALL = np.sort(glob.glob('*.c'))
  INC_ALL = np.sort(glob.glob('*.h'))
  SRC = ''
  OBJ = ''
  INC = ''
  for n in range(len(SRC_ALL)):
    SRC += '%s ' % SRC_ALL[n]
    OBJ += '%s.o ' % SRC_ALL[n][:-2]
  for n in range(len(INC_ALL)):
    INC += '%s ' % INC_ALL[n]

  # WRITE MAKEFILE
  mf = open('makefile', 'w')
  mf.write('CC = ' + host['COMPILER'] + '\n')
  mf.write('CCFLAGS = ' + C_FLAGS + ' ' + LIBRARIES + ' ' + INCLUDES + '\n')
  mf.write('LIB_FLAGS = ' + LIB_FLAGS + '\n')
  mf.write('CC_COMPILE = $(CC) $(CCFLAGS) -c' + '\n')
  mf.write('CC_LOAD = $(CC) $(CCFLAGS)' + '\n')
  mf.write('.c.o:' + '\n')
  mf.write('\t$(CC_COMPILE) $*.c' + '\n')
  mf.write('EXE = ' + EXECUTABLE + '\n')
  mf.write('all: $(EXE)' + '\n')
  mf.write('SRC = ' + SRC + '\n')
  mf.write('OBJ = ' + OBJ + '\n')
  mf.write('INC = ' + INC + '\n')
  mf.write('$(OBJ) : $(INC) makefile' + '\n')
  mf.write('$(EXE): $(OBJ) $(INC) makefile' + '\n')
  mf.write('\t$(CC_LOAD) $(OBJ) $(LIB_FLAGS) -o $(EXE)')
  mf.close()

  print("\n  COMPILING SOURCE\n")

  ncomp = 0
  first_error = 1
  if DEBUG:
    popen = subprocess.Popen(['make'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)
  else:
    popen = subprocess.Popen(['make','-j','10'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)

  for stdout_line in iter(popen.stdout.readline, ""):
    if stdout_line.rstrip()[-2:] == '.c':
      print("    [" + color.BOLD + color.BLUE +
            "%2d%%" % (100.*float(ncomp)/len(SRC_ALL)) + color.NORMAL +
            "] " + color.BOLD +
            stdout_line.rsplit('-c',1)[1].rstrip().lstrip().split('/')[-1] +
            color.NORMAL)
      ncomp += 1
  for stderr_line in iter(popen.stderr.readline, ""):
    if first_error == 1 and 'error' in stderr_line:
      print(color.WARNING + "\n  COMPILER ERROR\n" + color.NORMAL)
      first_error = 0
    print(stderr_line.rstrip())

  if first_error != 1:
    print(color.WARNING + "\n  COMPILATION FAILED\n" + color.NORMAL)
    sys.exit()

  obj_files = glob.glob('*.o')
  for f in obj_files:
    os.remove(f)
  os.remove('makefile')
  #os.rename(PATHS['SRC'] + 'bhlight', PATHS['BUILD'] + 'bhlight')

  print("\n  BUILD SUCCESSFUL")

  print("")

  sys.exit()

if len(sys.argv) != 2:
  print('ERROR: Format is')
  print('  python build.py [model]')
  sys.exit()
if (not os.path.isfile('model/' + sys.argv[1] + '.c')):
  print('ERROR Model %s does not exist' % sys.argv[1])
  sys.exit()

import shutil
shutil.copyfile('model/' + sys.argv[1] + '.c', 'model.c')
shutil.copyfile('model/' + sys.argv[1] + '.h', 'model.h')
build()


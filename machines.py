import os
import sys
machines = {}

def add_machine(name, compiler, c_flags, l_flags, gsl_dir):
  machine = {}
  machine['NAME'] = name
  machine['COMPILER'] = compiler
  machine['COMPILER_FLAGS'] = c_flags
  machine['LIB_FLAGS'] = l_flags
  machine['GSL_DIR'] = gsl_dir
  machines[name] = machine

def get_machine():
  for key in machines:
    if machines[key]['NAME'] in os.uname()[1]:
    #if os.uname()[1] == machines[key]['NAME']:
      return machines[key]
  print('UNKNOWN MACHINE ' + os.uname()[1])
  sys.exit()

add_machine(name='meade', 
            compiler='h5pcc',  
            c_flags='-O3 -Wall -Werror -fdiagnostics-color -fopenmp',
            l_flags='',
            gsl_dir='/home/brryan/Software/gsl')

add_machine(name='bh',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g -DVERSION=\"$(GIT_VERSION)\"',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='bh21',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g -DVERSION=\"$(GIT_VERSION)\"',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='bh27',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g -DVERSION=\"$(GIT_VERSION)\"',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='lmc',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g -DVERSION=\"$(GIT_VERSION)\"',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='stampede2',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g -DVERSION=\"$(GIT_VERSION)\"',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='/opt/apps/intel17/gsl/2.3')

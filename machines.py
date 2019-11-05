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

add_machine(name='elbert',
            compiler='h5pcc',
            c_flags='-Ofast -std=c99 -Wall -fopenmp',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

# IL Campus cluster
add_machine(name='ccc', # intel
            compiler='h5pcc',
            c_flags='-O3 -xCORE-AVX512 -std=c99 -Wall -qopenmp',
            l_flags='-lgsl -lgslcblas',
            gsl_dir='/usr/local/src/gsl/2.5')

add_machine(name='golubh',
            compiler='h5pcc',
            c_flags='-O3 -xCORE-AVX512 -std=c99 -Wall -qopenmp',
            l_flags='-lgsl -lgslcblas',
            gsl_dir='/usr/local/src/gsl/2.5')

add_machine(name='bh21',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='bh27',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g -D_DEFAULT_SOURCE',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='bh28',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='bh',
            compiler='/usr/lib64/mpi/gcc/openmpi/bin/h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g -D_DEFAULT_SOURCE',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='lmc',
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='')

add_machine(name='stampede2gcc', # gcc
            compiler='h5pcc',
            c_flags='-O3 -std=c99 -Wall -fopenmp -g',
            l_flags='-lm -lgsl -lgslcblas',
            gsl_dir='/opt/apps/gcc7_1/gsl/2.3')

add_machine(name='stampede2', # intel
            compiler='h5pcc',
            c_flags='-O3 -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -std=c99 -Wall -qopenmp -g',
            l_flags='-lgsl -lgslcblas',
            gsl_dir='/opt/apps/intel17/gsl/2.3')
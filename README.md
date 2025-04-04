# ipole
Polarized covariant radiative transfer in C, for imaging black hole accretion systems such as those being imaged by the Event Horizon Telescope. `ipole` is described in Moscibrodzka and Gammie [2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...43M/abstract).

This file details how to use `ipole` and includes a note about when it might not be be appropriate for your project.

# Prerequisites

The only prerequisites are the GNU Scientific Library (GSL) and HDF5. Specifically, the executable ```h5cc``` should be in your ```PATH```.

On Illinois EHT/BH cluster, this means loading the modules
```bash
$ module load gnu hdf5
```

On Linux machines, if you have root access, these can be installed to the system with
```bash
$ sudo apt install libhdf5-dev libgsl-dev
```
or
```bash
$ sudo dnf install hdf5-devel hdf5-static gsl-devel
```
The packages are straightforward to compile from source as well.

The native macOS version of ```clang``` does not support OpenMP parallelization, so on macOS with ```brew```, install
```bash
$ brew install llvm gsl hdf5
```
these should all be binary packages, no compiling should be required.

# Building

Then just build by running

```bash
$ make
```
with or without any of the usual ```make``` options, e.g., ```-j```.  You can 
make `ipole` into any destination directory by going there and then 
specifying the location of the makefile with ```-f```.  Note that
`ipole` does not support building in directories that contain spaces.

If you are using your own version of HDF5, you can specify where to find it with
```bash
$ make CC=/path/to/h5cc
```

... or compile with parallel HDF5:

```bash
$ make CC=/path/to/h5pcc
```

A particular fluid model or data format can be specified with
```MODEL=model_name```. For nearly all GRMHD fluid data, the ```iharm``` model
should be used -- this is the default if no model is specified.  After building
a model, be sure to run ```make clean``` before building a different one.

The optional argument ```NOTES=``` can be any string containing no spaces.
This string will be "baked-in" to the executable and printed out whenever the
program is run.

# Running

### Parameter files

To facilitate runs over many files with similar parameters,
`ipole` can read parameter files, specified with
```-par path/to/file.par```.

Parameter files should be lists of key value pairs, one per line.

```
# Example parameter
thetacam 90.0
```

Empty lines and lines beginning with the ```#``` character are ignored. A list of supported parameters for each model can be found in the ```example.par``` file in its directory. Expected keys that do not appear in the parameter file will be set to a backward-compatible default value (see, e.g., ```par.c``` for more information). Unexpected keys are silently ignored for forward compatibility.

### Command line

The same parameters can also be specified on the command line by using double-dash syntax,

```bash
$ ./ipole -par default.par --key1=value1 --key2=value2
$ ./ipole --freqcgs=230.e9 --MBH=6.2e9 --M_unit=1.e25 --thetacam=17 --dump=/path/to/dump.h5 --outfile=image.h5
```

In all cases, ```ipole``` reads parameters (and parameter files) in the order  they are specified, overwriting when encountering a repeated parameter name. Thus for most use cases the parameter file should be specified *first*.

There are two command-line specific flags:

```-unpol```: Only compute using approximate "unpolarized" transfer equation. This speeds up computation and is useful for (e.g.) fitting Munit.

```-quench```: Don't write any files to disk.

These options do not affect the rest of the arguments.

### Tracing

This version of ipole can output a file containing everything it used to compute
a particular pixel.  This includes the coordinate locations of each step along
the geodesic, as well as several scalars at each point (see the
[format doc](https://github.com/AFD-Illinois/docs/wiki/Trace-File-Output-Format)
for specifics).

To activate tracing for a pixel (i,j), measured from the bottom left with i
rightward along the x axis and j upward along y, run ipole with

```bash
$ ./ipole -par parameters.par --trace=1 --trace_i=i --trace_j=j
```

an alternative option ```trace_stride``` is available to trace every N pixels
across the entire image.

# Output

The ipole output format for images is documented
[here](https://github.com/AFD-Illinois/docs/wiki/Image-Format).
Converters are available to translate output into FITS images, or into the old
7-column text format (i,j,F,I,Q,U,V).

The ipole output format for traces is documented
[here](https://github.com/AFD-Illinois/docs/wiki/Trace-File-Output-Format)

# ipole limitations

```ipole``` is not a general-purpose imaging code. It has been designed, tested for,
and used for modeling of particular sources.  It may crash or
produce incorrect results in other use cases.  We strongly recommend you test
the code in an appropriate regime before you use it.

The ipole algorithm is described in Moscibrodzka and Gammie 
[2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...43M/abstract).  The
original version can be found [here](https://github.com/moscibrodzka/ipole).
Each radiative transfer step is based on an analytic solution to the
polarized transfer equation assuming constant coefficients. This
is intended to produce sensible results even when absorption or Faraday rotation
and conversion are large in a single step.  A comparison of
polarized relativistic radiative transport schemes can be found in
Prather [2023](https://ui.adsabs.harvard.edu/abs/2023ApJ...950...35P/abstract).

If you are imaging GRMHD simulations, you will find guidance on how a GRMHD
simulation and analysis pipeline fits together in the PATOKA pipeline paper,
Wong [2022](https://ui.adsabs.harvard.edu/abs/2022ApJS..259...64W/abstract).

```ipole``` treats synchrotron emission and absorption but not bremsstrahlung or Compton
scattering.  The transfer coefficients (emissivities, absorptivities, and rotativities)
are drawn from analytic fits presented
in Dexter [2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.462..115D/abstract) , Pandya 
[2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...822...34P/abstract), Pandya 
[2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...868...13P/abstract), and
Marszewski [2021](https://arxiv.org/abs/2108.10359), with a summary in the latter.  All 
coefficients assume that the electron distribution function is isotropic and that the frequency 
is large compared to the plasma frequency and cyclotron frequency. ipole implements fits for 
several electron distribution, including a thermal (Maxwell-Juttner) distribution, a power-law 
distribution, and a kappa distribution.

# Support and Contributing

We welcome contributions and pull requests from the community implementing new features or
issue fixes.  Not every new application fits the `ipole` scope, however, and we reserve
the right to choose which features to support.

You are welcome to file issues in this repository to ask questions of general interest and
report when `ipole` doesn't work for you. Please keep in mind that while bugs related to
`ipole`'s primary purpose will likely be given attention quickly, new uses and extensions
can only be supported on a best-effort basis.

ipole is an open source, community code that originated in the AFD group at Illinois
and in the Astrophysics group at Radboud.  Unfortunately we are not funded to provide support for all
users, so you may need to do some detective work yourself.  We will prioritize
fixing problems that are manifestly bugs, rather than problems encountered in
extending ipole into new parameter regimes.

# Funding

Funding for the development of ipole at Illinois was provided by the National
Science Foundation under multiple grants, including NSF AST 17-16327,
NSF OISE 17-43747, NSF AST 20-07936, and NSF AST 20-34306.


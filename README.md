# ipole
Polarized covariant radiative transfer in C, for imaging black hole accretion systems such as those being imaged by the Event Horizon Telescope. `ipole` is described in Moscibrodzka and Gammie [2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...43M/abstract).

This file details how to use ipole, along with a note about when it might (or might not) be appropriate for your project.

# Prerequisites

The only prerequisites are the GNU Scientific Library (GSL), and HDF5. Specifically,
the executable ```h5cc``` should be in your ```PATH```.

On Illinois BH cluster, this means loading the modules
```bash
$ module load gnu hdf5
```

On Linux machines, if you have root access these can be installed to the system with
```bash
$ sudo apt install libhdf5-dev libgsl-dev
```
or
```bash
$ sudo dnf install hdf5-devel hdf5-static gsl-devel
```
they are each straightforward to compile from source as well.

The native macOS version of ```clang``` does not support OpenMP parallelization.  So,
on macOS with ```brew```, install
```bash
$ brew install llvm gsl hdf5
```
these should all be binary packages, no compiling should be required.

# Building

Then just build by running

```bash
$ make
```
with or without any of the usual ```make``` options, e.g. ```-j```.  You can 
make ```ipole``` into any destination directory by going there, and then 
specifying the location of the makefile with ```-f```.  Note, though, that
```ipole``` does not support building in directories which contain spaces.

If you are using your own version of HDF5, you can specify where to find it with
```bash
$ make CC=/path/to/h5cc
```

Or compile with parallel HDF5:

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
```ipole``` can read parameter files, specified with
```-par path/to/file.par```.

Parameter files should be lists of key value pairs, one per line.

```
# Example parameter
thetacam 90.0
```

Empty lines and lines beginning with the ```#``` character are ignored.
A complete list of supported parameters for each model can be found in the
```example.par``` file in its directory. Expected keys
that do not appear in the parameter file will be set to a backward-compatible
default value (see ```par.c``` for more information). Unexpected keys are
silently ignored for forward-compatibility.

### Command line

The same parameters can also be specified on the command line by using
double-dash syntax,

```bash
$ ./ipole -par default.par --key1=value1 --key2=value2
$ ./ipole --freqcgs=230.e9 --MBH=6.2e9 --M_unit=1.e25 --thetacam=17 --dump=/path/to/dump.h5 --outfile=image.h5
```

In all cases, ```ipole``` reads parameters (and parameter files) in the order 
they are specified, overwriting when encountering a repeated parameter name.
Thus for most use cases the parameter file should be specified *first*.

There are two command-line specific flags:

```-unpol```: Only compute unpolarized transfer. This speeds up computation and
is useful for (e.g.) fitting an Munit.

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

# Before you use ipole for science

ipole is a powerful tool, but not a general-purpose imaging code. It is primarily designed to image RIAF accretion systems like M87* and SgrA*, specifically modeling only synchrotron emission from one of a number of electron energy distributions, without scattering. While it may work for your application, when used outside this domain (e.g. for very low temperatures, high optical depths, etc.) ipole may crash, or worse, produce arbitrarily incorrect results.

As a part of producing scientifically valid images, you should understand the domains of validity of the transfer coefficient fitting functions we use (Dexter [2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.462..115D/abstract) , Pandya et al. [2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...822...34P/abstract), Pandya et al. [2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...868...13P/abstract), Marszewski et al. 2021, submitted).  You should also have a passing understanding of the algorithm itself, described in (Moscibrodzka and Gammie [2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...43M/abstract)). If you are imaging GRMHD simulations, you should also understand the easy mistakes to make in imaging GRMHD results, and how the simulation and analysis pipeline fits together (e.g. Wong et al 2021, in prep).

ipole is provided to the community as-is, in the hope it will be useful.  We (the AFD Group at Illinois) are happy to answer questions and provide usage advice; however, development support in adapting the code for new use cases cannot be guaranteed.

# Contributing

We welcome contributions and pull requests from the community implementing new features or issue fixes, and we'll do our best to review and incorporate changes into the main version, to keep as many features and fixes available to as many people as possible.  Not every new application fits ipole's scope, however, and we reserve some discretion in choosing which features to incorporate (and consequently support).

Feel free to file issues in this repository where ipole doesn't work for your use case for any reason, or to ask questions which might be of general interest. Please keep in mind that while bugs will likely be given attention quickly, new uses and extensions must be supported on a best-effort basis.


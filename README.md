# ipole
Polarized covariant radiative transfer in C, for imaging black hole accretion systems such as those being imaged by the Event Horizon Telescope. `ipole` is described in Moscibrodzka and Gammie [2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...43M/abstract).

This file details how to use `ipole` and includes a note about when it might not be be appropriate for your project.

# Prerequisites

The only prerequisites are the GNU Scientific Library (GSL), and HDF5. Specifically, the executable ```h5cc``` should be in your ```PATH```.

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

# Before you use ipole for science

`ipole` is a powerful tool, but it is not a general-purpose imaging code. It is designed to image RIAF accretion systems like M87* and SgrA*, and it only models synchrotron emission and absorption using transfer coefficients from a particular electron distribution function. `ipole` does not treat scattering. While `ipole` may work for your use case, when used outside this domain (e.g., for very low temperatures, high optical depths, etc.) `ipole` may crash or produce misleading results.

As a part of producing scientifically valid images, you should understand the domains of validity of the transfer coefficient fitting functions that we use (Dexter [2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.462..115D/abstract), Pandya et al. [2016](https://ui.adsabs.harvard.edu/abs/2016ApJ...822...34P/abstract), Pandya et al. [2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...868...13P/abstract), Marszewski et al. 2021, submitted).  You should also have a passing understanding of the polarized transfer algorithm, which is described in (Moscibrodzka and Gammie [2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...43M/abstract)). If you are imaging GRMHD simulations, you should also understand the easy mistakes to make in imaging GRMHD results and how the simulation and analysis pipeline fits together (e.g., Wong et al 2021, in prep).

`ipole` is provided to the community as-is, in the hope it will be useful.  Although the AFD Group at Illinois, who maintain this respository, are happy to answer questions and provide usage advice, development support in adapting the code for new use cases cannot be guaranteed.

# Contributing

We welcome contributions and pull requests from the community implementing new features or issue fixes, and we will do our best to review and merge modifications so that the new features and fixes are available to as many people as possible.  Not every new application fits the `ipole` scope, however, and we reserve the right to choose which features to support.

You are welcome to file issues in this repository to ask questions of general interest and report when `ipole` doesn't work for you. Please keep in mind that while bugs related to the `ipole`'s primary purpose will likely be given attention quickly, new uses and extensions will only be supported on a best-effort basis.


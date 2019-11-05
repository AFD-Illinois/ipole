# ipole
Polarized covariant radiative transfer in C.

# Building

The only prerequisite is HDF5, specifically the executable ```h5pcc``` should 
be in your ```PATH```. Then just build by running

```bash
$ make
```
A particular fluid model or data format can be specified with
```MODEL=model_name```. For all fluid data produced from Illinois codes after
September 2018, the default ```iharm``` model should be used.

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
```

In all cases, ```ipole``` reads parameters (and parameter files) in the order 
they are specified, overwriting when encountering a repeated parameter name.
Thus for most use cases the parameter file should be specified *first*.

There are two command-line specific flags:

```-unpol```: Only compute unpolarized transfer. This speeds up computation and
is useful for (e.g.) fitting an Munit.

```-quench```: Don't write any files to disk.

These options do not affect the rest of the arguments.

# ipole
polarized covariant radiative transport code

# building

You can build the ```ipole``` executable by using the python build script as follows
```bash
$ python build.py [MODEL] [NOTES]
```
where ```[MODEL]``` specifies the fluid data format. For all fluid data produced from Illinois codes after September 2018, the ```iharm``` model should be used. The optional argument ```[NOTES]``` can be any string (containing no spaces). This string will be "baked-in" to the executable and printed out whenever the program is run. 

# running

```ipole``` can be run using either parameter files or command line arguments. ```ipole``` also accepts the following flags:

```-unpol``` Only compute unpolarized transfer. This speeds up computation and is useful for (e.g.) fitting an Munit.

```-quench``` Don't write any files to disk.

These options may be provided at any time in any order and do not affect the rest of the arguments.

### command line arguments and parameter files

To facilitate runs over many files with similar parameters, ```ipole``` can read parameter files, specified by consecutive the command line arguments ```-par path/to/par.file```.

Parameter files should be lists of key value pairs, written one line at a time. Empty lines and lines beginning with the ```#``` character are ignored. A full list of accepted parameters can be found by reading the ```par.c``` file. An example parameter file can be found at [the bottom](#example-parameter-file) of this document. Expected keys that do not appear in the parameter file will be set to default value (see ```par.c``` for more information). Unexpected keys that appear in the parameter file will be silently ignored.

Values in a parameter file can be overridden on the command line by using double-dash syntax as 
```bash
$ ./ipole -par default.par --key1=value1 --key2=value2
$ ./ipole --key1=value1 -par default.par --key2=value2
```

In all cases, ```ipole``` reads parameters in the order they are specified, overwriting when encountering a repeated parameter name.


## program loading order

Knowing this order may be useful for adding new variables that should be loaded / overwritten before / during / after other parameters are set.

1. ```par.c:load_par(...)``` may be called if triggered by the command line arguments ```-par path/to/par/file```. This function populates the ```Params``` struct with various user-defined variables (often similar to/a superset of values you would pass at runtime) such as inclination angle, M_unit, and so on. If this function is called, the ```params``` data structure is marked as having been loaded.
2. ```model.c:parse_input(...)``` is called. This function is responsible for populating model-defined variables such as M_unit, freqcgs, and so on. These variables tend to live within the model.c/h files (as static globals) and thus should not be referenced from without. If the params data structure is marked as loaded from above, then this function ignores command line input and sets all model variables according to the params structure. (TODO in the future, it might be nice to allow command line arguments to override parameter file values.)
3. ```model.c:init_model(...)``` is called. This function is responsible for setting up the initial data structures associated with the fluid dump file(s) and grid geometry. By the time this function has returned, all initial dump files should be loaded into memory. Before returning, this function may call ```load_*_data(...)``` below.
4. ```model.c:load_*_data(...)``` may be called on initial load or during the ```update_data``` stage, if ipole is running in the slow-light configuration.


## examples

### example parameter file
```
# example.par  --  an example parameter file
# Note that lines beginning with the '#' character and
# empty lines are ignored by the parser. 

dump /path/to/fluid/dump.h5
thetacam 90
phicam 0 
freqcgs 230e.9
MBH 1.e8
M_unit 1.e40
outfile data/image.h5


```




# ipole
polarized covariant radiative transport code

# building

You can build the ```ipole``` executable by using the python build script as follows
```bash
$ python build.py [MODEL]
```
where ```[MODEL]``` specifies the fluid data format. For all fluid data produced from Illinois codes after September 2018, the ```iharm``` model should be used.

# running

ipole can be run using either parameter files or command line arguments. 

### command line arguments

If no parameter file is specified, then ipole reads arguments in the following format

```bash
./ipole thetacam freqcgs Mbh M_unit path/to/dump counterjet

# as example
./ipole 17 230.e9 1.e8 2.3e21 /data/bh28-fs1/dumps/iharm/mad/0.9375_256x128x128/dumps/dump_00000000.h5 0
```

```thetacam``` should be given in degrees away from the polar axis.

```freqcgs``` is the desired frequency for the image (given in Hz).

```Mbh``` is the mass of the central black hole in units of Msun (overwritten if the dump comes from a GRRMHD run)

```M_unit``` is the mass unit for the simulation in cgs (overwritten if the dump compes from a GRRMHD run)

```path/to/dump``` is the absolute or relative path to the fluid dump file

```counterjet``` is one of {0,1,2} and specifies if certain parts of the domain should be "switched off" for emission. If counterjet == 1, only emission from X[2] > 0.5 is allowed. If counterjet == 2, only emission from X[2] < 0.5 is allowed. Otherwise, emission from all parts of the domain is allowed.

### parameter files

TODO


## program loading order

Knowing this order may be useful for adding new variables that should be loaded / overwritten before / during / after other parameters are set.

1. ```par.c:load_par(...)``` may be called if triggered by the command line arguments ```-par path/to/par/file```. This function populates the ```Params``` struct with various user-defined variables (often similar to/a superset of values you would pass at runtime) such as inclination angle, M_unit, and so on. If this function is called, the ```params``` data structure is marked as having been loaded.
2. ```model.c:parse_input(...)``` is called. This function is responsible for populating model-defined variables such as M_unit, freqcgs, and so on. These variables tend to live within the model.c/h files (as static globals) and thus should not be referenced from without. If the params data structure is marked as loaded from above, then this function ignores command line input and sets all model variables according to the params structure. (TODO in the future, it might be nice to allow command line arguments to override parameter file values.)
3. ```model.c:init_model(...)``` is called. This function is responsible for setting up the initial data structures associated with the fluid dump file(s) and grid geometry. By the time this function has returned, all initial dump files should be loaded into memory. Before returning, this function may call ```load_*_data(...)``` below.
4. ```model.c:load_*_data(...)``` may be called on initial load or during the ```update_data``` stage, if ipole is running in the slow-light configuration.



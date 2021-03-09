# Crustal Deformation and Fault Mechanics Group, Stanford University
# 
## cdfm_fdra2c
CDFM magic happens here...

## To Compile:
This code does not (yet) have a proper `.configure` or even `autoconf` setup just yet. It probably requires a bit of refactoring, to move some functionality from code to compiler. In particular, the `util/` directory should probably be moved up a level and treated as a separate package. This is easy enough, excet that codes/headers in this sub-package are referenced by relative directory in the code, ie `include util/include/my_header.hpp` (or something like that) so that the compiler flag is `-I external`, and `util` is located in the `external/` subdirectory. 

As suggested, the `util/` code should nominally be elevated level with `fdra2c`, and then facilitate a simple `include my_headder.hpp` via compiler flags.


For now, this should work (maybe with some minor changes) to compile on Mazama HPC:

```
module purge;module load gnu/8 openmpi_3/;make clean;HMMVP_DIR=external/hmmvp UTIL_DIR=external/util make mode=p opt=-O3 -j 4
```

The nested `util/` package compiled, but it will probalby need to be recompiled for a given platform, and minor adjustment to`Makefile` might be necessary, though I hope that the next step will be to configure for `autotools`, `cmake`, or another builder (`pymake`?).

### NOTE:
The main significance that this compile is for Mazama HPC is the choice of compiler, specifically in this case, `gnu/8.3.x` and `openmpi/3.x`.
  

# OptimizeMPI.jl

This repository contains source code for the following poster at Julia Con 2020 : 
"Calling a parallel simulation code from Julia". The aim is to demonstrating to call
a distributed memory parallel program from julia for doing parameter optimization. I want
to share my small expertiences about that

# Directories
- `squared_norm` : a very basic example using the squared L₂ norm
- `heat` : control problem on the 2D heat equation using finite differencies
- `examples` : small code snippets used in the notebook
# Requirements 

- CMake 
- MPI library
- a Fortran 2003 Compiler (tested against Gfortran 10.1.0)

# Compilation
- go into `<directory>` (`heat` or `squared_norm`)
```bash
cd heat
```
- If your CMake is recent enough run the following
```bash
cmake -B build && for i in test, install; do make -j 4 -C build $i ; done
```
It must install a *shared library* in the current directory (.so, .dylib or .dll)
Note : If your cmake does not suport -B, you could use the legacy way

```bash
mkdir build && cd build && cmake .. && make && make test && make install && cd ..
```

# Running examples

```bash
mpirun -np 9 julia scr_optim.jl
```
# Poster
[![Click_here](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/aitzkora/OptimizeMPI.jl/blob/master/calling_a_parallel_code.ipynb?flush_cache=true)

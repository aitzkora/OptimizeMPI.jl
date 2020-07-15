# OptimizeMPI.jl

This repository contains source code for a poster at Julia Con 2020 :
"Calling a parallel simulation code from Julia". The aim is to present a basic tutorial
explaining how to call a distributed memory parallel program from Julia and optimize
some of its parameters.

# directories

- `squared_norm` : a very basic example using the squared $L_2$ norm
- `heat` : control problem on the 2D heat equation using finite differencies
- `examples` : small code snippets used in the notebook

# requirements 

- CMake 
- MPI library
- a Fortran 2003 Compiler (tested against Gfortran 10.1.0)
- and gcc for snippets
- Julia (tested with 1.5.0-rc1, but should work for versions >=1.x) with the 
following packages installed : 
 - `MPI.jl` : necessary
 - `MPIClusterManagers.jl` : optional, just to be more comfortable when running examples.

# compilation
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
# running examples
If you have installed the package `MPIClusterManagers.jl` on your Julia version, you could run

```bash
julia cluster_optim.jl
```
otherwise, `mpirun` must launch julia for yourself

```bash
mpirun -np 9 julia scr_optim.jl
```
For clarity, a sequential version of the heat equation is implemented in Julia, in the file
`heat_seq.jl`

# poster
[![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/aitzkora/OptimizeMPI.jl/blob/master/calling_a_parallel_code.ipynb?flush_cache=true)
# Troubleshooting

Sometimes, the compilation of `MPI.jl` or `MPIClusterManagers.jl` does not choose the MPI library that
you want. Remember you can set some variable in the `startup.jl` file like
```julia
ENV["JULIA_MPI_C_LIBRARIES"] = "-L/usr/lib/openmpi/ -lmpi"
ENV["JULIA_MPI_Fortran_INCLUDE_PATH"] = "-I/usr/include"
ENV["JULIA_MPI_PATH"] = "/usr/"
```

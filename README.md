# OptimizeMPI.jl

This repository contains source code for the following poster at Julia Con 2020 : 
"Calling a parallel simulation code from Julia"

# Directories
- `squared_norm` : very basic example using the squared L₂ norm
- `heat` : control problem on the 2D heat equation using finite differencies

# Requirements 

- cmake 
- MPI library
- a fortran compiler

# Compilation
- go into <directory> (`heat` or `squared_norm`)
```bash
cd heat
```
- run cmake 
```bash
mkdir build && cd build cmake && make && make test
```

# Running examples

```bash
julia scr_optim.jl
```

*todo*: add the poster and documentation using a ipython notebook

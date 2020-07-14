import MPI
include("test_optim.jl")
res, x_min = find_optim(100);

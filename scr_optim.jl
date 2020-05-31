using MPIClusterManagers, Distributed
manager = MPIManager(np=4)
addprocs(manager)
println("Added procs $(procs())")
@everywhere import MPI


println("running test_optim")
@mpi_do manager begin 
    include("test_optim.jl");
    res, x_min = find_optim(100);
end

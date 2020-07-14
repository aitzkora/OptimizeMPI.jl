using MPIClusterManagers, Distributed
manager = MPIManager(np=4)
addprocs(manager)
println("Added procs $(procs())")
@everywhere import MPI

@mpi_do manager begin 
    include("test_optim.jl");
    find_optim(100,false);
end

using MPIClusterManagers, Distributed
manager = MPIManager(np=4)
addprocs(manager)
println("Added procs $(procs())")
@everywhere import MPI

@mpi_do manager begin 
    include("heat_par.jl");
    find_optim(10, false);
end
